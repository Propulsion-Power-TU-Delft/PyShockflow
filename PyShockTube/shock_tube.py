import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
import csv
from PyShockTube.riemann_problem import RiemannProblem
from PyShockTube.roe_scheme import RoeScheme_Base, RoeScheme_Generalized
from PyShockTube.muscl_hancock import MusclHancock
from PyShockTube.fluid import FluidIdeal, FluidReal
from PyShockTube.euler_functions import *



class ShockTube:
    def __init__(self, config):
        """
        Initializes the problem with space and time arrays, along with additional fluid properties.

        Parameters
        ----------
        config: configuration file object

        Returns
        -------
        None
        """
        print("\n" + "=" * 80)
        print(" " * 25 + "🚀  WELCOME TO PYSHOCKTUBE 🚀")
        print(" " * 18 + "Fluid Dynamics Simulation for Shock Tubes")
        print("=" * 80)
        
        print()  
        print("=" * 80)
        print(" "*32 + "SIMULATION DATA")
        # print("Length of the domain [m]:                    %.2e" % (x[-1] - x[0]))
        # print("Number of points:                            %i" %len(x))
        # print("Final time instant [s]:                      %.2e" % (t[-1]))
        # print("Fluid name:                                  %s" % fluid_props[0])
        # print("Fluid treatment:                             %s" % fluid_props[1])
        
        self.config = config
        self.topology = self.config.getTopology()
        
        # geometry
        self.nNodes = self.config.getNumberOfPoints()
        self.length = self.config.getLength()
        self.xNodes = np.linspace(0, self.length, self.nNodes)
        self.dx = self.xNodes[1]-self.xNodes[0]
        
        # fluid initial states
        self.pressureLeft = self.config.getPressureLeft()
        self.pressureRight = self.config.getPressureRight()
        self.densityLeft = self.config.getDensityLeft()
        self.densityRight = self.config.getDensityRight()
        self.velocityLeft = self.config.getVelocityLeft()
        self.velocityRight = self.config.getVelocityRight()
        
        # compute dt
        Smax = np.sqrt(1.4*np.max([self.pressureLeft, self.pressureRight])/np.min([self.densityLeft, self.densityRight]))+np.max([self.velocityLeft, self.velocityRight])  # brutal approximation max eigenvalue
        CFLmax = self.config.getCFLMax() 
        dtMax = CFLmax * self.dx / Smax
        self.timeMax = self.config.getTimeMax()
        nt = int(self.timeMax/dtMax)        
        self.timeVec = np.linspace(0, self.timeMax, nt)
        self.dt = self.timeVec[1]-self.timeVec[0]
        self.nTime = len(self.timeVec)
        self.fluid_name = self.config.getFluidName()
        self.fluid_model = self.config.getFluidModel()
        if self.fluid_model.lower()=='ideal':
            self.gmma = self.config.getFluidGamma()
            self.Rgas = self.config.getGasRConstant()
            print("Fluid cp/cv ratio [-]:                       %.2e" %self.gmma)
            print("Fluid gas constant [J/kgK]:                  %.2e" %self.Rgas)
            self.fluid = FluidIdeal(self.gmma,self.Rgas)
        elif self.fluid_model.lower()=='real':
            self.fluid = FluidReal(self.fluid_name)
        
        self.BCtype = self.config.getBoundaryConditions()
            
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
        
        if self.topology.lower()=='default':
            self.areaTube = np.zeros_like(self.xNodesVirt)+1
        elif self.topology.lower()=='nozzle':
            self.areaTube = self.readNozzleFile(self.xNodesVirt, self.config.getNozzleFilePath())
            print('The simulation topology selected is nozzle. Reading the coordinates')
        else:
            raise ValueError('Unknown topology type')
        
        self.dAreaTude_dx = np.gradient(self.areaTube, self.xNodesVirt)
        
        


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
        

    def InitialConditionsLeftRight(self):
        """
        Initialize the conditions based on initial state, defined by right and left values
        """
        initialConditions = {'Density': np.array([self.densityLeft, self.densityRight]), 'Velocity': np.array([self.velocityLeft, self.velocityRight]), 'Pressure': np.array([self.pressureLeft, self.pressureRight])}
        
        print("Initial L/R density values [kg/m3]:          (%.2e, %.2e)" %(initialConditions['Density'][0], initialConditions['Density'][1]))
        print("Initial L/R velocity values [m/s]:           (%.2e, %.2e)" %(initialConditions['Velocity'][0], initialConditions['Velocity'][1]))
        print("Initial L/R pressure values [bar]:           (%.2e, %.2e)" %(initialConditions['Pressure'][0]/1e5, initialConditions['Pressure'][1]/1e5))

        initialConditions['Energy'] = self.fluid.ComputeStaticEnergy_p_rho(initialConditions['Pressure'], initialConditions['Density'])
        for name in self.solutionNames:
            self.solution[name][:, 0] = self.CopyInitialState(initialConditions[name][0], initialConditions[name][1])
    

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
        xInterface = self.config.getInterfaceLocation()
        assert(xInterface>0), f"The interface must be located within 0 and {self.length}"
        assert(xInterface<self.length), f"The interface must be located within 0 and {self.length}"
        
        f = np.zeros_like(self.xNodesVirt)
        for i in range(len(self.xNodesVirt)):
            if self.xNodesVirt[i] <= xInterface:
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


    def PlotNodeSolution(self, iNode, folder_name = None, file_name = None):
        """
        Plot the solution at a specified node location over time
        """
        fig, ax = plt.subplots(2, 2, figsize=(12, 8))
        ax[0, 0].plot(self.timeVec[:], self.solution['Density'][iNode, :], '-C0o', ms=2)
        ax[0, 0].set_ylabel(r'Density')

        ax[0, 1].plot(self.timeVec[:], self.solution['Velocity'][iNode, :], '-C1o', ms=2)
        ax[0, 1].set_ylabel(r'Velocity')

        ax[1, 0].plot(self.timeVec[:], self.solution['Pressure'][iNode, :], '-C2o', ms=2)
        ax[1, 0].set_ylabel(r'Pressure')

        ax[1, 1].plot(self.timeVec[:], self.solution['Energy'][iNode, :], '-C3o', ms=2)
        ax[1, 1].set_ylabel(r'Energy')

        fig.suptitle('Location %.3f' % self.xNodes[iNode])

        for row in ax:
            for col in row:
                col.set_xlabel('t')
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


    def SetBoundaryConditions(self, it):
        """
        Set the correct boundary condition type (`reflective`, `transparent`, or `periodic`)
        """
        

        if it==0:
            print("Boundary Conditions Left:                    %s" %self.BCtype[0])
            print("Boundary Conditions Right:                   %s" %self.BCtype[1])
            print("="*80)

        if self.BCtype[0].lower()=='reflective':
            self.SetReflectiveBoundaryConditions(it, 'left')
        elif self.BCtype[0].lower()=='transparent':
            self.SetTransparentBoundaryConditions(it, 'left')
        elif self.BCtype[0].lower()=='periodic':
            self.SetPeriodicBoundaryConditions(it, 'left')
        else:
            raise ValueError("Unknown boundary condition type on the left")
        
        if self.BCtype[1].lower()=='reflective':
            self.SetReflectiveBoundaryConditions(it, 'right')
        elif self.BCtype[1].lower()=='transparent':
            self.SetTransparentBoundaryConditions(it, 'right')
        elif self.BCtype[1].lower()=='periodic':
            self.SetPeriodicBoundaryConditions(it, 'right')
        else:
            raise ValueError("Unknown boundary condition type on the right")
        
        # update also the conservative variable arrays based on what has been done on the primitive
        self.solutionCons['u1'][:, it], self.solutionCons['u2'][:, it], self.solutionCons['u3'][:, it] = (
                    GetConservativesFromPrimitives(self.solution['Density'][:, it], self.solution['Velocity'][:, it], self.solution['Pressure'][:, it], self.fluid))


    def SetReflectiveBoundaryConditions(self, iTime, location):
        """
        Set reflective BC at time `iTime`
        """
        if location=='left':
            self.solution['Density'][0, iTime] = self.solution['Density'][1, iTime]
            self.solution['Velocity'][0, iTime] = -self.solution['Velocity'][1, iTime]
            self.solution['Pressure'][0, iTime] = self.solution['Pressure'][1, iTime]
            self.solution['Energy'][0, iTime] = self.solution['Energy'][1, iTime]
        elif location=='right':
            self.solution['Density'][-1, iTime] = self.solution['Density'][-2, iTime]
            self.solution['Velocity'][-1, iTime] = -self.solution['Velocity'][-2, iTime]
            self.solution['Pressure'][-1, iTime] = self.solution['Pressure'][-2, iTime]
            self.solution['Energy'][-1, iTime] = self.solution['Energy'][-2, iTime]
        else:
            raise ValueError('Unknown location specified')
            
    

    def SetTransparentBoundaryConditions(self, iTime, location):
        """
        Set transparent BC at time `iTime`
        """
        if location=='left':
            self.solution['Density'][0, iTime] = self.solution['Density'][1, iTime]
            self.solution['Velocity'][0, iTime] = self.solution['Velocity'][1, iTime]
            self.solution['Pressure'][0, iTime] = self.solution['Pressure'][1, iTime]
            self.solution['Energy'][0, iTime] = self.solution['Energy'][1, iTime]
        elif location=='right':
            self.solution['Density'][-1, iTime] = self.solution['Density'][-2, iTime]
            self.solution['Velocity'][-1, iTime] = self.solution['Velocity'][-2, iTime]
            self.solution['Pressure'][-1, iTime] = self.solution['Pressure'][-2, iTime]
            self.solution['Energy'][-1, iTime] = self.solution['Energy'][-2, iTime]
        else:
            raise ValueError('Unknown location specified')
        
        
        
    def SetPeriodicBoundaryConditions(self, iTime, location):
        """
        Set periodic BC at time `iTime`
        """
        if location=='left':
            self.solution['Density'][0, iTime] = self.solution['Density'][-2, iTime]
            self.solution['Velocity'][0, iTime] = self.solution['Velocity'][-2, iTime]
            self.solution['Pressure'][0, iTime] = self.solution['Pressure'][-2, iTime]
            self.solution['Energy'][0, iTime] = self.solution['Energy'][-2, iTime]
        elif location=='right':
            self.solution['Density'][-1, iTime] = self.solution['Density'][1, iTime]
            self.solution['Velocity'][-1, iTime] = self.solution['Velocity'][1, iTime]
            self.solution['Pressure'][-1, iTime] = self.solution['Pressure'][1, iTime]
            self.solution['Energy'][-1, iTime] = self.solution['Energy'][1, iTime]
        else:
            raise ValueError('Unknown location specified')


    def SolveSystem(self, high_order=False, limiter='Van Albada'):
        """
        Solve the equations explicitly in time (forward Euler) using a certain flux_method (`Godunov`, `Roe`, `WAF`). high_order
        specifies if applying or not high order reconstruction with limiters. At the moment only type one is working -> simply
        impose high_order=True
        """
        flux_method = self.config.getNumericalScheme()
        print()
        print("="*80)
        print(" "*33 + "START SOLVER")
        print("Numerical flux method: %s" %(flux_method))
        print("MUSCL reconstruction: %s" %high_order)
        # if high_order:
        #     print("Limiter: %s" %(limiter))
        print()

        # short aliases
        cons = self.solutionCons
        prim = self.solution
        dx, dt = self.dx, self.dt
        
        # time-steps loop
        for it in range(1, self.nTime):
            print('Time step: %i of %i' %(it, self.nTime))

            # compute fluxes on every internal interface
            flux = np.zeros((self.nNodes+1, 3))
            for iFace in range(flux.shape[0]):
                flux[iFace, :] = self.ComputeFluxVector(iFace, iFace+1, it-1, flux_method, high_order, limiter)
            
            if self.topology.lower()=='nozzle':
                source = self.ComputeSourceTerms(it-1)
            else:
                source = np.zeros((self.nNodesHalo,3))
            
            # update the conservatives for every element
            for iNode in range(1, self.nNodesHalo-1):
                cons['u1'][iNode, it] = cons['u1'][iNode, it-1] + dt/dx * ((flux[iNode-1, 0] - flux[iNode, 0]) + source[iNode, 0]*dx)
                cons['u2'][iNode, it] = cons['u2'][iNode, it-1] + dt/dx * ((flux[iNode-1, 1] - flux[iNode, 1]) + source[iNode, 1]*dx)
                cons['u3'][iNode, it] = cons['u3'][iNode, it-1] + dt/dx * ((flux[iNode-1, 2] - flux[iNode, 2]) + source[iNode, 2]*dx)

            # update the primitives
            prim['Density'][1:-1, it], prim['Velocity'][1:-1, it], prim['Pressure'][1:-1, it], prim['Energy'][1:-1, it] = \
                GetPrimitivesFromConservatives(cons['u1'][1:-1, it], cons['u2'][1:-1, it], cons['u3'][1:-1, it], self.fluid)
            
            # set boundary conditions to update the ghost points for the new iteration
            self.SetBoundaryConditions(it)

        print(" "*34 + "END SOLVER")
        print("="*80)
    
    
    def ComputeSourceTerms(self, it):
        source_terms = np.zeros((self.nNodesHalo,3))
        source_terms[:,0] = - self.solution['Density'][:, it] * self.solution['Velocity'][:, it]*self.dAreaTude_dx/self.areaTube
        source_terms[:,1] = - (self.solution['Density'][:, it] * self.solution['Velocity'][:, it]**2)*self.dAreaTude_dx/self.areaTube
        source_terms[:,2] = - self.solution['Velocity'][:, it] *(self.solution['Energy'][:, it] + self.solution['Pressure'][:, it])*self.dAreaTude_dx/self.areaTube
        return source_terms
        


    def ComputeFluxVector(self, il, ir, it, flux_method, high_order, limiter):
        """
        Compute the flux vector at the interface between grid points `il` and `ir`, using a certain `flux_method`.
        """
        
        # flow reconstruction if high_order=True
        if (high_order and il>2 and ir<self.nNodesHalo-2):
            rhoL, uL, pL, rhoR, uR, pR = self.MUSCL_VanAlbada_Reconstruction(il, ir, it)  # currently working
            # rhoL, uL, pL, rhoR, uR, pR = self.MUSCL_Reconstruction(il, ir, it, limiter)  # at the moment creates oscillations with every limiter for some reasons
        else:
            rhoL = self.solution['Density'][il, it]
            rhoR = self.solution['Density'][ir, it]
            uL = self.solution['Velocity'][il, it]
            uR = self.solution['Velocity'][ir, it]
            pL = self.solution['Pressure'][il, it]
            pR = self.solution['Pressure'][ir, it]            
        
        # flux calculation
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


    def MUSCL_VanAlbada_Reconstruction(self, il, ir, it):
        """
        MUSCL approach with Van Albada limiter. Formulation taken from pag. 110 of "Computational Fluid Dynamics book, by Blazek", where kappa=0,
        and the reconstruction with limiter is embedded in a single formula
        """
        # states left, left minus 1, right, right plus one
        U_l = np.array([self.solution['Density'][il, it], self.solution['Velocity'][il, it], self.solution['Pressure'][il, it]])
        U_lm = np.array([self.solution['Density'][il-1, it], self.solution['Velocity'][il-1, it], self.solution['Pressure'][il-1, it]])
        U_r = np.array([self.solution['Density'][ir, it], self.solution['Velocity'][ir, it], self.solution['Pressure'][ir, it]])
        U_rp = np.array([self.solution['Density'][ir+1, it], self.solution['Velocity'][ir+1, it], self.solution['Pressure'][ir+1, it]])

        # unlimited jumps
        aR = U_rp-U_r
        bR = U_r-U_l
        aL = U_r-U_l
        bL = U_l-U_lm

        def func(a, b, eps=1e-6):
            y = (a*(b**2+eps)+b*(a**2+eps)) / (a**2 + b**2 +2*eps) 
            return y
        
        # limited slopes
        deltaR = func(aR, bR)
        deltaL = func(aL, bL)

        # reconstruct left and right states
        U_l_rec = U_l+0.5*deltaL
        U_r_rec = U_r-0.5*deltaR

        return U_l_rec[0], U_l_rec[1], U_l_rec[2], U_r_rec[0], U_r_rec[1], U_r_rec[2]


    def ShowAnimation(self):
        """
        Show animation of the results at all time instants
        """
        ni, nt = self.solution['Density'].shape
        mach = np.zeros((ni, nt))
        for i in range(ni):
            for t in range(nt):
                mach[i,t] = self.solution['Velocity'][i,t]/self.fluid.ComputeSoundSpeed_p_rho(self.solution['Pressure'][i,t], self.solution['Density'][i,t])
        def plot_limits(f, extension=0.05):
            max = f.max()
            min = f.min()
            left = min-(max-min)*extension
            right = max+(max-min)*extension
            return left, right
        
        if self.config.showAnimation():
            fig, ax = plt.subplots(2, 2, figsize=(12, 8))
            density_limits = plot_limits(self.solution['Density'])
            velocity_limits = plot_limits(self.solution['Velocity'])
            pressure_limits = plot_limits(self.solution['Pressure'])
            energy_limits = plot_limits(self.solution['Energy'])
            mach_limitis = plot_limits(mach)
            for it in range(self.nTime):
                for row in ax:
                    for col in row:
                        col.cla()
                ax[0, 0].plot(self.xNodesVirt, self.solution['Density'][:, it], '-C0o', ms=2)
                ax[0, 0].set_ylabel(r'Density [kg/m3]')
                ax[0, 0].set_ylim(density_limits)

                ax[0, 1].plot(self.xNodesVirt, self.solution['Velocity'][:, it], '-C1o', ms=2)
                ax[0, 1].set_ylabel(r'Velocity [m/s]')
                ax[0, 1].set_ylim(velocity_limits)

                ax[1, 0].plot(self.xNodesVirt, self.solution['Pressure'][:, it], '-C2o', ms=2)
                ax[1, 0].set_ylabel(r'Pressure [Pa]')
                ax[1, 0].set_ylim(pressure_limits)

                # ax[1, 1].plot(self.xNodesVirt, self.solution['Energy'][:, it], '-C3o', ms=2)
                # ax[1, 1].set_ylabel(r'Energy')
                # ax[1, 1].set_ylim(energy_limits)
                ax[1, 1].plot(self.xNodesVirt, mach[:, it], '-C3o', ms=2)
                ax[1, 1].set_ylabel(r'Mach [-]')
                ax[1, 1].set_ylim(mach_limitis)

                fig.suptitle('Time %.3e [s]' % self.timeVec[it])

                for row in ax:
                    for col in row:
                        col.set_xlabel('x')
                        col.grid(alpha=.3)
                plt.pause(1e-3)
            
            
            massflow = self.solution['Density'][-1,:]*self.solution['Velocity'][-1,:]
            plt.figure()
            plt.plot(self.timeVec, massflow)
            plt.xlabel('Time [s]')
            plt.ylabel('Mass flow exit [kg/s]')
            plt.grid(alpha=.3)
            
            id_throat = np.where(self.areaTube==np.min(self.areaTube))[0][0]
            
            plt.figure()
            plt.plot(self.timeVec, mach[-1,:], label='exit')
            plt.plot(self.timeVec, mach[id_throat,:], label='throat')
            plt.xlabel('Time [s]')
            plt.ylabel('Mach [-]')
            plt.grid(alpha=.3)
            plt.legend()
            
            ni = mach.shape[0]
            plt.figure()
            plt.plot(self.timeVec, self.solution['Pressure'][2*ni//4,:], label='L/2')
            plt.plot(self.timeVec, self.solution['Pressure'][3*ni//4,:], label='3L/4')
            plt.xlabel('Time [s]')
            plt.ylabel('Pressure [Pa]')
            plt.grid(alpha=.3)
            plt.legend()
            plt.show()
    

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
    

    def SaveSolution(self):
        """
        Save the full object as a pickle for later use
        """
        folder_name = self.config.getOutputFolder()
        os.makedirs(folder_name, exist_ok=True)
        file_name = self.config.getOutputFileName()
        full_path = folder_name+'/'+file_name+'_%.6f.pik' %(self.timeMax)
        with open(full_path, 'wb') as file:
            pickle.dump(self, file)
        print('Pickle object with full solution saved to ' + full_path + ' !')


    def SaveNodeSolutionToCSV(self, iNode, timeInstants, folder_name, file_name):
        """
        Save the 'Pressure' and 'Temperature' array from the solution to a CSV file.
        """
        file_path = folder_name + '/' + file_name + '.dat'
        pressure_data = self.solution['Pressure'][iNode, :]  # Extract the pressure data (1D array)
        density_data = self.solution['Density'][iNode, :]  # Extract the density data (1D array)
        temperature_data = self.fluid.ComputeTemperature_p_rho(pressure_data, density_data)

        with open(file_path, 'w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            for value in range(len(timeInstants)):
                writer.writerow([timeInstants[value], pressure_data[value], temperature_data[value]])
            # writer.writerow(pressure_data)  # Write each row of the array to the CSV file
            # writer.writerow(timeInstants)  # Write each row of the array to the CSV file

        print(f"Pressure and temperature data saved to {file_path}!")

    def MUSCL_Reconstruction(self, il, ir, it, limiter):
        """
        MUSCL reconstruction with limiter. Currently not working, since it creates oscillations with every limiter
        """

        # states left, left minus 1, right, right plus one
        U_l = np.array([self.solution['Density'][il, it], self.solution['Velocity'][il, it], self.solution['Pressure'][il, it]])
        U_lm = np.array([self.solution['Density'][il-1, it], self.solution['Velocity'][il-1, it], self.solution['Pressure'][il-1, it]])
        U_r = np.array([self.solution['Density'][ir, it], self.solution['Velocity'][ir, it], self.solution['Pressure'][ir, it]])
        U_rp = np.array([self.solution['Density'][ir+1, it], self.solution['Velocity'][ir+1, it], self.solution['Pressure'][ir+1, it]])
        
        def compute_jump_ratio(num, den):
            r = np.zeros_like(num)
            for i in range(len(num)):
                if np.abs(num[i])<=1e-8:
                    num[i] = 0.
                    den[i] = 1.
                elif (num[i]>=1e-8 and np.abs(den[i])<1e-8):
                    num[i] = 1.
                    den[i] = 1.
                elif (num[i]<-1e-8 and np.abs(den[i])<1e-8):
                    num[i] = -1.
                    den[i] = 1.
            return num/den
        
        num = (U_r-U_l)
        den = (U_l-U_lm)
        r_l = compute_jump_ratio(num, den)

        num = (U_rp-U_r)
        den = (U_r-U_l)
        r_r = compute_jump_ratio(num, den)

        Phi_l = self.Compute_Limiter(r_l, limiter)
        Phi_r = self.Compute_Limiter(r_r, limiter)

        U_l_rec = U_l + 0.5*Phi_l*(U_r-U_l)
        U_r_rec = U_r - 0.5*Phi_r*(U_r-U_l)

        return U_l_rec[0], U_l_rec[1], U_l_rec[2], U_r_rec[0], U_r_rec[1], U_r_rec[2]


    def Compute_Limiter(self, r_vec, limiter):
        """
        Compute the limiter values. Currently not working.
        """
        psi = np.zeros(3)
        for i in range(len(r_vec)):
            r = r_vec[i]

            if limiter.lower() == 'van albada':
                psi[i] = (r**2+r)/(1+r**2)

            elif limiter.lower() == 'van leer':
                psi[i] = (r+np.abs(r))/(1+np.abs(r))

            elif limiter.lower() == 'min mod':
                psi[i] = np.maximum(0, np.minimum(1, r))

            elif limiter.lower() == 'superbee':
                psi[i] = np.maximum(0, np.minimum(2*r, 1), np.minimum(r, 2))

            elif limiter.lower() == 'none':
                psi[i] = 1 
            else:
                raise ValueError('Limiter not recognized!')
        
        return psi
    
    def readNozzleFile(self, xTube, filepath):
        nozzleData = np.loadtxt(filepath, skiprows=1, delimiter=',', dtype=float)
        nozzleX = nozzleData[:,0]
        nozzleArea = nozzleData[:,1]
        
        # Linear interpolation with external filling set to 1
        interpolatedNozzleArea = np.interp(xTube, nozzleX, nozzleArea, left=1, right=1)
        return interpolatedNozzleArea
        























