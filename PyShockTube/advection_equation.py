import numpy as np
import matplotlib.pyplot as plt
import os
import pickle


class AdvectionEquation:
    def __init__(self, x, t, cfl):
        """
        Initialize the problem with space and time arrays.
        """
        self.xNodes = x
        self.nNodes = len(x)
        self.dx = x[1]-x[0]
        self.timeVec = t
        self.dt = t[1]-t[0]
        self.nTime = len(t)
        self.nInterfaces = self.nNodes+1
        self.cfl = cfl
        self.a_speed = self.dx*self.cfl/self.dt
        self.gmma = 1.4
        self.GenerateVirtualGeometry()
        self.u = np.zeros((self.nNodesHalo, self.nTime))


    def GenerateVirtualGeometry(self):
        """
        Halo-nodes at the boundaries
        """
        self.nNodesHalo = self.nNodes+2
        self.xNodesVirt = np.zeros(self.nNodesHalo)
        self.xNodesVirt[1:-1] = self.xNodes
        self.xNodesVirt[0] = self.xNodes[0]-self.dx
        self.xNodesVirt[-1] = self.xNodes[-1]+self.dx
        self.blockage = np.ones_like(self.xNodesVirt)
        self.cell_volumes = np.zeros(self.nNodesHalo)+self.dx
        self.int_surface = np.zeros(self.nInterfaces)+1

    def InstantiateInitialCondition(self, wave_type, wave_length):
        """
        Initialize the initial condition, on the first wave_length of x-domain
        """
        L = self.xNodes[-1]-self.xNodes[0]
        if wave_type=='square':
            for i in range(len(self.u)):
                if self.xNodesVirt[i] >0 and self.xNodesVirt[i] < wave_length:
                    self.u[i, 0] = 1
        elif wave_type=='sine':
            for i in range(len(self.u)):
                if self.xNodesVirt[i] >0 and self.xNodesVirt[i] < wave_length:
                    nwaves = 1
                    self.u[i, 0] = np.sin(nwaves*np.pi*self.xNodes[i]/wave_length)
        elif wave_type=='triangular':
            for i in range(len(self.u)):
                if self.xNodesVirt[i] >0 and self.xNodesVirt[i] < wave_length/2:
                    slope = 1/(wave_length/2)
                    self.u[i, 0] = self.xNodes[i]*slope
                elif self.xNodesVirt[i] >=wave_length/2 and self.xNodesVirt[i] < wave_length:
                    slope = -1/(wave_length/2)
                    self.u[i, 0] = 1+(self.xNodes[i]-wave_length/2)*slope
        elif wave_type=='saw-tooth':
            for i in range(len(self.u)):
                if self.xNodesVirt[i] >0 and self.xNodesVirt[i] <= wave_length:
                    slope = 1/(wave_length)
                    self.u[i, 0] = self.xNodes[i]*slope
        else:
            raise ValueError('unknown wave_type parameter')
        
        self.u_analytic = self.u.copy()
        

    def PlotSolution(self, iTime, folder_name = None, file_name = None):
        """
        Plot the solution at time instant element iTime
        :param iTime: element index in time array
        :param folder_name:
        :param file_name:
        :return:
        """
        plt.figure()
        plt.plot(self.xNodesVirt, self.u_analytic[:, iTime], label='Reference')
        try: plt.plot(self.xNodesVirt, self.u[:, iTime], label='%s' %self.method)
        except: pass

        plt.title('Time %.3f' % self.timeVec[iTime])
        plt.xlabel('x')
        plt.ylabel('u')
        plt.grid(alpha=.3)

        if file_name is not None and folder_name is not None:
            os.makedirs(folder_name, exist_ok=True)
            plt.savefig(folder_name + '/' + file_name + '.pdf', bbox_inches='tight')
    

    def SetBoundaryConditions(self, BCs, it):
        """
        Set the correct boundary condition type. Only periodic for advection equation
        """
        self.BCtype = BCs
        if BCs=='periodic':
            self.SetPeriodicBoundaryConditions(it)
        else:
            raise ValueError("Unknown boundary condition type")
        


    def ComputeBlockagedGeometry(self, b):
        
        #shrink the cell volumes
        for iNode in range(1, self.nNodesHalo-1):
            iBlock = iNode-1
            self.cell_volumes[iNode] *= b[iBlock]

        # shrink the surfaces area
        for iInt in range(self.nInterfaces):
            if iInt>0 and iInt<self.nInterfaces-1:  # no need to treat the extremes, since the blockage is not there probably
                block = 0.5*(b[iInt-1]+b[iInt])
                self.int_surface[iInt] *= block
    
    
    def SetPeriodicBoundaryConditions(self, iTime):
        """
        Set the boundary conditions, making use of the halo nodes (index 0 and -1 along axis 0), and then update the conservative variables
        :parameter iTime: time instant index of the simulation
        """
        self.u[0, iTime] = self.u[-2, iTime]
        self.u[-1, iTime] = self.u[1, iTime]

    def SolveSystem(self, method):
        """
        Solve the advection using a certain method. Choose between: Lax-Wendroff, Godunov-upwind, Godunov-centred,
        FORCE, Lax-Friedrichs
        """
        self.method = method
        dx, dt = self.dx, self.dt
        alpha = self.GetAlphaScheme(method)
        for it in range(1, self.nTime):
            print('Time step: %i of %i' %(it, self.nTime))
            for ix in range(1, self.nNodesHalo-1):
                if ix==1:
                    fL = 0.5*(1+2*alpha*self.cfl)*(self.a_speed*self.u[ix-1, it-1]) + 0.5*(1-2*alpha*self.cfl)*(self.a_speed*self.u[ix, it-1])
                else:
                    fL = fR  # use of conservation 
                fR = 0.5*(1+2*alpha*self.cfl)*(self.a_speed*self.u[ix, it-1]) + 0.5*(1-2*alpha*self.cfl)*(self.a_speed*self.u[ix+1, it-1])



                flux_net = fL-fR
                self.u[ix, it] = self.u[ix, it-1] + dt/dx*flux_net
            self.SetBoundaryConditions(self.BCtype, it)
            self.u_analytic[:, it] = self.SolveAnalytical(it)
    
    def SolveAnalytical(self, it):
        """
        For every time step, also find the analytical solution advecting the initial wave
        """
        u_pad = np.concatenate((self.u[:, 0], self.u[:, 0]))
        u_sol = self.u[:, 0].copy()
        for i in range(self.nNodesHalo):
            time = self.timeVec[it]
            deltaX = time*self.a_speed
            shifted_x_idx = int(deltaX/self.dx)
            u_sol[i] = u_pad[len(u_sol)+i-shifted_x_idx]
        return u_sol 
    
    def GetAlphaScheme(self, method):
        if method=='Lax-Wendroff':
            return 0.5
        elif method=='Godunov-upwind':
            return 1/2/self.cfl
        elif method=='Godunov-centred':
            return 1
        elif method=='FORCE':
            return 1/4/self.cfl**2 *(1+self.cfl**2)
        elif method=='Lax-Friedrichs':
            return 1/2/self.cfl**2


    def ShowAnimation(self):
        """
        Show animation of the results for all time instants
        """
        plt.figure()
        for it in range(self.nTime):
            plt.cla()
            plt.plot(self.xNodesVirt, self.u_analytic[:, it], label='Reference')
            plt.plot(self.xNodesVirt, self.u[:, it], label='%s' %self.method)
            plt.xlabel('x')
            plt.ylabel('u')
            plt.grid(alpha=0.3)
            plt.legend()

            plt.title('Time %.3f' % self.timeVec[it])
            plt.pause(1e-3)
    

    def SaveSolution(self, folder_name, file_name):
        """
        Save the full object
        """
        os.makedirs(folder_name, exist_ok=True)
        full_path = folder_name+'/'+file_name+'.pik'
        with open(full_path, 'wb') as file:
            pickle.dump(self, file)
        print('Solution saved to ' + full_path + ' !')























