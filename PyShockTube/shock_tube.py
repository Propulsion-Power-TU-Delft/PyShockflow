import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
from NumericalCodes.riemann_problem import RiemannProblem
from NumericalCodes.roe_scheme import RoeScheme_Base, RoeScheme_Generalized
from NumericalCodes.muscl_hancock import MusclHancock
from NumericalCodes.fluid import FluidIdeal, FluidReal
from NumericalCodes.euler_functions import *



class ShockTube:
    def __init__(self, x, t, *fluid_props):
        """
        Initializes the problem with space and time arrays, along with additional fluid properties.

        Parameters
        ----------
        x: array_like, x coordinates of the 1D-grid

        t: array_like, time instants of the simulation

        *fluid_props: property of the fluid to use (fluid name, fluid model, fluid gamma)

        Returns
        -------
        None
        """
        self.xNodes = x
        self.nNodes = len(x)
        self.dx = x[1]-x[0]
        self.timeVec = t
        self.dt = t[1]-t[0]
        self.nTime = len(t)
        self.fluid_name = fluid_props[0]
        if fluid_props[1].lower()=='ideal':
            assert(len(fluid_props)>=3)
            self.fluid_model = 'ideal'
            self.gmma = fluid_props[2]
            self.fluid = FluidIdeal(self.gmma)
        elif fluid_props[1].lower()=='real':
            assert(len(fluid_props)>=2)
            self.fluid_model = 'real'
            self.fluid = FluidReal(self.fluid_name)
        self.GenerateVirtualGeometry()


    def GenerateVirtualGeometry(self):
        """
        Generate the virtual geometry consisting of halo nodes for boundary conditions

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.nNodesHalo = self.nNodes+2
        self.xNodesVirt = np.zeros(self.nNodesHalo)
        self.xNodesVirt[1:-1] = self.xNodes
        self.xNodesVirt[0] = self.xNodes[0]-self.dx
        self.xNodesVirt[-1] = self.xNodes[-1]+self.dx


    def InstantiateSolutionArrays(self):
        """
        Instantiate the containers for the solutions. The first dimension is space, the second is time.
        """
        self.solutionNames = ['Density', 'Velocity', 'Pressure', 'Energy']
        self.solution = {}
        for name in self.solutionNames:
            self.solution[name] = np.zeros((self.nNodesHalo, self.nTime))


    def InstantiateSolutionArraysConservatives(self):
        """
        Instantiate the containers for the solutions. The first dimension is space, the second is time.
        """
        self.solutionConsNames = ['u1', 'u2', 'u3']
        self.solutionCons = {}
        for name in self.solutionConsNames:
            self.solutionCons[name] = np.zeros((self.nNodesHalo, self.nTime))
        

    def InitialConditionsLeftRight(self, dictIn):
        """
        Initialize the conditions based on initial state, defined by right and left values
        """
        dictIn['Energy'] = self.fluid.ComputeStaticEnergy_p_rho(dictIn['Pressure'], dictIn['Density'])
        for name in self.solutionNames:
            self.solution[name][:, 0] = self.CopyInitialState(dictIn[name][0], dictIn[name][1])
    

    def InitialConditionsArrays(self, dictIn):
        """
        Initialize the conditions based on initial state, defined by arrays
        """
        dictIn['Energy'] = dictIn['Pressure'] / (self.gmma - 1) / dictIn['Density']
        for name in self.solutionNames:
            self.solution[name][:, 1:-1] = dictIn[name]


    def CopyInitialState(self, fL, fR):
        """
        Given left and right values, copy these values along the x-axis
        :param fL:
        :param fR:
        :return:
        """
        xmean = 0.5 * (self.xNodesVirt[-1] + self.xNodesVirt[0])
        f = np.zeros_like(self.xNodesVirt)
        for i in range(len(self.xNodesVirt)):
            if self.xNodesVirt[i] <= xmean:
                f[i] = fL
            else:
                f[i] = fR
        return f


    def PlotSolution(self, iTime, folder_name = None, file_name = None):
        """
        Plot the solution at time instant element iTime
        """
        fig, ax = plt.subplots(2, 2, figsize=(12, 8))
        ax[0, 0].plot(self.xNodesVirt, self.solution['Density'][:, iTime], '-C0o', ms=2)
        ax[0, 0].set_ylabel(r'Density')

        ax[0, 1].plot(self.xNodesVirt, self.solution['Velocity'][:, iTime], '-C1o', ms=2)
        ax[0, 1].set_ylabel(r'Velocity')

        ax[1, 0].plot(self.xNodesVirt, self.solution['Pressure'][:, iTime], '-C2o', ms=2)
        ax[1, 0].set_ylabel(r'Pressure')

        ax[1, 1].plot(self.xNodesVirt, self.solution['Energy'][:, iTime], '-C3o', ms=2)
        ax[1, 1].set_ylabel(r'Energy')

        fig.suptitle('Time %.3f' % self.timeVec[iTime])

        for row in ax:
            for col in row:
                col.set_xlabel('x')
                col.grid(alpha=.3)

        if file_name is not None and folder_name is not None:
            os.makedirs(folder_name, exist_ok=True)
            plt.savefig(folder_name + '/' + file_name + '.pdf', bbox_inches='tight')
    

    def PlotConservativeSolution(self, iTime, folder_name = None, file_name = None):
        """
        Plot the conservative variables at time instant element iTime
        """
        fig, ax = plt.subplots(1, 3, figsize=(15, 4))
        ax[0].plot(self.xNodesVirt, self.solutionCons['u1'][:, iTime], '-C0o', ms=2)
        ax[0].set_ylabel(r'$\rho$')

        ax[1].plot(self.xNodesVirt, self.solutionCons['u2'][:, iTime], '-C1o', ms=2)
        ax[1].set_ylabel(r'$\rho u$')

        ax[2].plot(self.xNodesVirt, self.solutionCons['u3'][:, iTime], '-C2o', ms=2)
        ax[2].set_ylabel(r'$\rho e$')

        fig.suptitle('Time %.3f' % self.timeVec[iTime])

        for col in ax:
            col.set_xlabel('x')
            col.grid(alpha=.3)

        if file_name is not None and folder_name is not None:
            os.makedirs(folder_name, exist_ok=True)
            plt.savefig(folder_name + '/' + file_name + '_conservatives.pdf', bbox_inches='tight')


    def SetBoundaryConditions(self, BCs, it):
        """
        Set the correct boundary condition type (`reflective`, `transparent`, or `periodic`)
        """
        self.BCtype = BCs
        if BCs.lower()=='reflective':
            self.SetReflectiveBoundaryConditions(it)
        elif BCs.lower()=='transparent':
            self.SetTransparentBoundaryConditions(it)
        elif BCs.lower()=='periodic':
            self.SetPeriodicBoundaryConditions(it)
        else:
            raise ValueError("Unknown boundary condition type")
        
        # update also the conservative variable arrays based on what has been done on the primitive
        self.solutionCons['u1'][:, it], self.solutionCons['u2'][:, it], self.solutionCons['u3'][:, it] = (
                    GetConservativesFromPrimitives(self.solution['Density'][:, it], self.solution['Velocity'][:, it], self.solution['Pressure'][:, it], self.fluid))


    def SetReflectiveBoundaryConditions(self, iTime):
        """
        Set reflective BC at time `iTime`
        """
        self.solution['Density'][0, iTime] = self.solution['Density'][1, iTime]
        self.solution['Density'][-1, iTime] = self.solution['Density'][-2, iTime]

        self.solution['Velocity'][0, iTime] = -self.solution['Velocity'][1, iTime]
        self.solution['Velocity'][-1, iTime] = -self.solution['Velocity'][-2, iTime]

        self.solution['Pressure'][0, iTime] = self.solution['Pressure'][1, iTime]
        self.solution['Pressure'][-1, iTime] = self.solution['Pressure'][-2, iTime]

        self.solution['Energy'][0, iTime] = self.solution['Energy'][1, iTime]
        self.solution['Energy'][-1, iTime] = self.solution['Energy'][-2, iTime]
    

    def SetTransparentBoundaryConditions(self, iTime):
        """
        Set transparent BC at time `iTime`
        """
        self.solution['Density'][0, iTime] = self.solution['Density'][1, iTime]
        self.solution['Density'][-1, iTime] = self.solution['Density'][-2, iTime]

        self.solution['Velocity'][0, iTime] = self.solution['Velocity'][1, iTime]
        self.solution['Velocity'][-1, iTime] = self.solution['Velocity'][-2, iTime]

        self.solution['Pressure'][0, iTime] = self.solution['Pressure'][1, iTime]
        self.solution['Pressure'][-1, iTime] = self.solution['Pressure'][-2, iTime]

        self.solution['Energy'][0, iTime] = self.solution['Energy'][1, iTime]
        self.solution['Energy'][-1, iTime] = self.solution['Energy'][-2, iTime]
    

    def SetPeriodicBoundaryConditions(self, iTime):
        """
        Set periodic BC at time `iTime`
        """
        self.solution['Density'][0, iTime] = self.solution['Density'][-2, iTime]
        self.solution['Density'][-1, iTime] = self.solution['Density'][1, iTime]

        self.solution['Velocity'][0, iTime] = self.solution['Velocity'][-2, iTime]
        self.solution['Velocity'][-1, iTime] = self.solution['Velocity'][1, iTime]

        self.solution['Pressure'][0, iTime] = self.solution['Pressure'][-2, iTime]
        self.solution['Pressure'][-1, iTime] = self.solution['Pressure'][1, iTime]

        self.solution['Energy'][0, iTime] = self.solution['Energy'][-2, iTime]
        self.solution['Energy'][-1, iTime] = self.solution['Energy'][1, iTime]


    def SolveSystem(self, flux_method):
        """
        Solve the equations explicitly in time using a certain flux_method (`Godunov`, `Roe`, `WAF`)
        """
        cons = self.solutionCons
        prim = self.solution
        dx, dt = self.dx, self.dt
        for it in range(1, self.nTime):
            print('Time step: %i of %i' %(it, self.nTime))
            for ix in range(1, self.nNodesHalo-1):
                if ix==1: # only for the first node compute left and right fluxes
                    fluxVec_left = self.ComputeFluxVector(ix-1, ix, it-1, flux_method)
                    fluxVec_right = self.ComputeFluxVector(ix, ix+1, it-1, flux_method)
                else:
                    fluxVec_left = fluxVec_right  # make use of the previously calculated flux (conservative approach)
                    fluxVec_right = self.ComputeFluxVector(ix, ix+1, it-1, flux_method)

                fluxVec_net = fluxVec_left-fluxVec_right
                
                cons['u1'][ix, it] = cons['u1'][ix, it-1] + dt/dx*fluxVec_net[0]
                cons['u2'][ix, it] = cons['u2'][ix, it-1] + dt/dx*fluxVec_net[1]
                cons['u3'][ix, it] = cons['u3'][ix, it-1] + dt/dx*fluxVec_net[2]

                prim['Density'][ix, it], prim['Velocity'][ix, it], prim['Pressure'][ix, it], prim['Energy'][ix, it] = \
                    GetPrimitivesFromConservatives(cons['u1'][ix, it], cons['u2'][ix, it], cons['u3'][ix, it], self.fluid)
                
            self.SetBoundaryConditions(self.BCtype, it)


    def ComputeFluxVector(self, il, ir, it, flux_method):
        """
        Compute the flux vector at the interface between grid points `il` and `ir`, using a certain `flux_method`
        """
        rhoL = self.solution['Density'][il, it]
        rhoR = self.solution['Density'][ir, it]
        uL = self.solution['Velocity'][il, it]
        uR = self.solution['Velocity'][ir, it]
        pL = self.solution['Pressure'][il, it]
        pR = self.solution['Pressure'][ir, it]
        
        if flux_method.lower()=='godunov':
            if self.fluid_model!='ideal':
                raise ValueError('Godunov scheme is available only for ideal gas model')
            nx, nt = 51, 51
            x = np.linspace(-self.dx/2, self.dx/2, nx)
            t = np.linspace(0, self.dt, nt)
            riem = RiemannProblem(x, t)
            riem.InitializeState([rhoL, rhoR, uL, uR, pL, pR])
            riem.InitializeSolutionArrays()
            riem.ComputeStarRegion()
            riem.Solve(space_domain='interface', time_domain='global') # compute Riemann solution only at x=0, but on all time instants
            rho, u, p = riem.GetSolutionInTime()
            u1, u2, u3 = GetConservativesFromPrimitives(rho, u, p, self.fluid)
            u1AVG, u2AVG, u3AVG = np.sum(u1)/len(u1), np.sum(u2)/len(u2), np.sum(u3)/len(u3)
            flux = EulerFluxFromConservatives(u1AVG, u2AVG, u3AVG, self.fluid) 
        
        elif flux_method.lower=='waf':
            if self.fluid_model!='ideal':
                raise ValueError('WAF scheme is available only for ideal gas model')
            nx, nt = 51, 51
            x = np.linspace(-self.dx/2, self.dx/2, nx)
            t = np.linspace(0, self.dt, nt)
            riem = RiemannProblem(x, t)
            riem.InitializeState([rhoL, rhoR, uL, uR, pL, pR])
            riem.InitializeSolutionArrays()
            riem.ComputeStarRegion()
            riem.Solve(space_domain='global', time_domain='interface') # compute Riemann solution only at deltaT/2, but on whole space-domain
            rho, u, p = riem.GetSolutionInSpace()
            u1, u2, u3 = GetConservativesFromPrimitives(rho, u, p, self.fluid)
            u1AVG, u2AVG, u3AVG = np.sum(u1)/len(u1), np.sum(u2)/len(u2), np.sum(u3)/len(u3)
            flux = EulerFluxFromConservatives(u1AVG, u2AVG, u3AVG, self.fluid)
        
        elif flux_method.lower()=='roe':
            if self.fluid_model=='ideal':
                roe = RoeScheme_Base(rhoL, rhoR, uL, uR, pL, pR, self.fluid)
                roe.ComputeAveragedVariables()
                roe.ComputeAveragedEigenvalues()
                roe.ComputeAveragedEigenvectors()
                roe.ComputeWaveStrengths()
                flux = roe.ComputeFlux()
            elif self.fluid_model=='real':
                roe = RoeScheme_Generalized(rhoL, rhoR, uL, uR, pL, pR, self.fluid)
                roe.ComputeAveragedVariables()
                roe.ComputeAveragedEigenvalues()
                roe.ComputeWaveStrengths()
                flux = roe.ComputeFlux()
            else:
                raise ValueError('Unknown fluid model. Select ideal or real')
        
        elif flux_method.lower()=='muscl-hancock':
            if self.fluid_model!='ideal':
                raise ValueError('MUSCL-Hancock scheme is available only for ideal gas model')
            # be careful when you are at the border, since in reality you would need two halo nodes, not just one
            if il>1 and ir<self.nNodesHalo-2:
                rhoLL, rhoRR = self.solution['Density'][il-1, it], self.solution['Density'][ir+1, it]
                uLL, uRR = self.solution['Velocity'][il-1, it], self.solution['Velocity'][ir+1, it]
                pLL, pRR = self.solution['Pressure'][il-1, it], self.solution['Pressure'][ir+1, it]
            else: # no extrapolation for extreme cells over the halo nodes
                rhoLL, rhoRR = rhoL, rhoR
                uLL, uRR = uL, uR
                pLL, pRR = pL, pR
            mhck = MusclHancock(rhoLL, rhoL, rhoR, rhoRR, uLL, uL, uR, uRR, pLL, pL, pR, pRR, self.dx)
            mhck.ReconstructInterfaceValues()
            mhck.EvolveInterfaceValues(self.dx, self.dt)
            flux = mhck.ComputeRoeFlux()
        
        else:
            raise ValueError('Unknown flux method')
        
        return flux

    def ShowAnimation(self):
        """
        Show animation of the results at all time instants
        """
        fig, ax = plt.subplots(2, 2, figsize=(12, 8))
        for it in range(self.nTime):
            for row in ax:
                for col in row:
                    col.cla()
            ax[0, 0].plot(self.xNodesVirt, self.solution['Density'][:, it], '-C0o', ms=2)
            ax[0, 0].set_ylabel(r'Density')

            ax[0, 1].plot(self.xNodesVirt, self.solution['Velocity'][:, it], '-C1o', ms=2)
            ax[0, 1].set_ylabel(r'Velocity')

            ax[1, 0].plot(self.xNodesVirt, self.solution['Pressure'][:, it], '-C2o', ms=2)
            ax[1, 0].set_ylabel(r'Pressure')

            ax[1, 1].plot(self.xNodesVirt, self.solution['Energy'][:, it], '-C3o', ms=2)
            ax[1, 1].set_ylabel(r'Energy')

            fig.suptitle('Time %.3f' % self.timeVec[it])

            for row in ax:
                for col in row:
                    col.set_xlabel('x')
                    col.grid(alpha=.3)
            plt.pause(1e-3)
    
    def ShowContourAnimation(self):
        """
        Show contour animation of the results for all time instants
        """
        xgrid = np.reshape(self.xNodes, (len(self.xNodes), 1))
        xgrid = np.hstack((xgrid, xgrid))
        ygrid = np.reshape(np.zeros_like(self.xNodes), (len(self.xNodes), 1))
        length = self.xNodes[-1]-self.xNodes[0]
        AR = 4
        height = length/AR
        ygrid = np.hstack((ygrid, ygrid+height))

        def convert2D(f):
            fgrid = np.reshape(f, (len(self.xNodes), 1))
            fgrid = np.hstack((fgrid, fgrid))
            return fgrid

        fig, ax = plt.subplots(2, 2, figsize=(8, 5))
        for it in range(self.nTime):
            for row in ax:
                for col in row:
                    col.cla()
            rhoGrid = convert2D(self.solution['Density'][1:-1, it])
            uGrid = convert2D(self.solution['Velocity'][1:-1, it])
            pGrid = convert2D(self.solution['Pressure'][1:-1, it])
            eGrid = convert2D(self.solution['Energy'][1:-1, it])
           
            ax[0, 0].contourf(xgrid, ygrid, rhoGrid, cmap='jet', levels=20)
            ax[0, 0].set_title(r'Density')

            ax[0, 1].contourf(xgrid, ygrid, uGrid, cmap='jet', levels=20)
            ax[0, 1].set_title(r'Velocity')

            ax[1, 0].contourf(xgrid, ygrid, pGrid, cmap='jet', levels=20)
            ax[1, 0].set_title(r'Pressure')

            ax[1, 1].contourf(xgrid, ygrid, eGrid, cmap='jet', levels=20)
            ax[1, 1].set_title(r'Energy')

            fig.suptitle('Time %.3f' % self.timeVec[it])

            for row in ax:
                for col in row:
                    col.set_aspect('equal')
            plt.pause(1e-3)
    

    def SaveSolution(self, folder_name, file_name):
        """
        Save the full object as a pickle for later use
        """
        os.makedirs(folder_name, exist_ok=True)
        full_path = folder_name+'/'+file_name+'.pik'
        with open(full_path, 'wb') as file:
            pickle.dump(self, file)
        print('Solution save to ' + full_path + ' !')























