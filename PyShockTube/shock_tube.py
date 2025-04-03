import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
import csv
import sys
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
        self.config = config
        self.topology = self.config.getTopology()
        self.fluid_name = self.config.getFluidName()
        self.fluid_model = self.config.getFluidModel()
        
        if self.fluid_model.lower()=='ideal':
            self.gmma = self.config.getFluidGamma()
            self.Rgas = self.config.getGasRConstant()
            self.fluid = FluidIdeal(self.gmma,self.Rgas)
        elif self.fluid_model.lower()=='real':
            self.fluid = FluidReal(self.fluid_name)
        
        # fluid initial states
        self.pressureLeft = self.config.getPressureLeft()
        self.pressureRight = self.config.getPressureRight()
        try:
            self.densityLeft = self.config.getDensityLeft()
            self.densityRight = self.config.getDensityRight()
        except:
            temperatureLeft = self.config.getTemperatureLeft()
            temperatureRight = self.config.getTemperatureRight()
            self.densityLeft = self.fluid.ComputeDensity_p_T(self.pressureLeft, temperatureLeft)
            self.densityRight = self.fluid.ComputeDensity_p_T(self.pressureRight, temperatureRight)
        self.velocityLeft = self.config.getVelocityLeft()
        self.velocityRight = self.config.getVelocityRight()
        
        # geometry
        self.length = self.config.getLength()
        self.nNodes = self.config.getNumberOfPoints()
        xNodes = self.GeneratePhysicalGeometry(self.length, self.nNodes)
        self.nNodes = len(xNodes)
        self.GenerateVirtualGeometry(xNodes)
        self.prepareOutputPaths()
        
        # Time related information
        self.cflMax = self.config.getCFLMax()
        self.timeMax = self.config.getTimeMax()
        
        # Boundary Conditions
        self.BCtype = self.config.getBoundaryConditions()    
        print("Boundary Conditions Left:                    %s" %self.BCtype[0])
        print("Boundary Conditions Right:                   %s" %self.BCtype[1])
        print("="*80)
        
        # Print info
        print("\n" + "=" * 80)
        print(" " * 25 + "🚀  WELCOME TO PYSHOCKTUBE 🚀")
        print(" " * 18 + "Fluid Dynamics Simulation for Shock Tubes")
        print("=" * 80)
        print()  
        print("=" * 80)
        print(" "*32 + "SIMULATION DATA")
        print("Length of the domain [m]:                    %.2e" % self.length)
        print("Number of points:                            %i" % self.nNodes)
        print("Final time instant [s]:                      %.2e" % self.timeMax)
        print("Fluid name:                                  %s" % self.fluid_name)
        print("Fluid treatment:                             %s" % self.fluid_model)
        if self.fluid_model.lower()=='ideal':
            print("Fluid cp/cv ratio [-]:                       %.2e" %self.gmma)
            print("Fluid gas constant [J/kgK]:                  %.2e" %self.Rgas)
    
    
    def prepareOutputPaths(self):
        self.outputFolder = self.config.getOutputFolder()
        os.makedirs(self.outputFolder, exist_ok=True)
        
        self.subFolder = self.outputFolder + '/' + self.config.getOutputFileName() + '_NX_%i' % (self.nNodes)
        dum = self.subFolder
        counter = 1
        while os.path.exists(dum):
            dum = f"{self.subFolder}_{counter}"
            counter += 1
        
        self.subFolder = dum
        os.makedirs(self.subFolder, exist_ok=True)
            
        
                

    def GeneratePhysicalGeometry(self, length, nodes):
        meshRefined = self.config.isMeshRefined()
        if meshRefined is False:
            xNodes = np.linspace(0, length, nodes)
        else:
            refinementCoords = self.config.getRefinementBoundaries()
            print("Mesh is refined between the two boundaries [m]: ", refinementCoords)
            
            pointsRefinement = self.config.getNumberPointsRefinement()
            totalPoints = nodes
            pointsOutside = totalPoints-pointsRefinement
            lengthUpstream = refinementCoords[0]
            lengthDownstream = length-refinementCoords[1]
            pointsUpstream = int(pointsOutside*(lengthUpstream)/(lengthUpstream+lengthDownstream))
            pointsDownstream = pointsOutside-pointsUpstream
            
            # case in which the refinement is internal
            if pointsDownstream>0 and pointsUpstream>0:
                xUpstream = np.linspace(0, refinementCoords[0], pointsUpstream+1)
                xRefinement = np.linspace(refinementCoords[0], refinementCoords[1], pointsRefinement+1)
                xDownstream = np.linspace(refinementCoords[1], length, pointsDownstream)
                if self.config.adaptMeshRefinementExtremities():
                    xUpstream = self.ComputeStretchedGridPoints(xUpstream, xRefinement, 'upstream')
                    xDownstream = self.ComputeStretchedGridPoints(xDownstream, xRefinement, 'downstream')
                xNodes = np.concatenate((xUpstream[0:-1], xRefinement[0:-1], xDownstream))
                
            elif pointsUpstream>0 and pointsDownstream==0: # the refinement finish with the end of the domain
                xUpstream = np.linspace(0, refinementCoords[0], pointsUpstream+1)
                xRefinement = np.linspace(refinementCoords[0], refinementCoords[1], pointsRefinement+1)
                if self.config.adaptMeshRefinementExtremities():
                    xUpstream = self.ComputeStretchedGridPoints(xUpstream, xRefinement, 'upstream')  
                xNodes = np.concatenate((xUpstream[0:-1], xRefinement))
            
            elif pointsUpstream==0 and pointsDownstream>0: # the refinement starts with the domain
                xRefinement = np.linspace(refinementCoords[0], refinementCoords[1], pointsRefinement+1)
                xDownstream = np.linspace(refinementCoords[1], length, pointsDownstream)
                if self.config.adaptMeshRefinementExtremities():
                    xDownstream = self.ComputeStretchedGridPoints(xDownstream, xRefinement, 'downstream')
                xNodes = np.concatenate((xRefinement[0:-1], xDownstream))
            
            else:
                raise ValueError('The refinement is ill-positioned. Please locate it internally to the domain, or at one of the extremities')
                
        return xNodes  
    
    
    def ComputeGridSpacing(self, xNodes):
        dx = np.zeros_like(xNodes)
        dx[0] = xNodes[1]-xNodes[0]
        for i in range(1,len(dx)-1):
            dx[i] = (xNodes[i+1]-xNodes[i])/2 + (xNodes[i]-xNodes[i-1])/2
        dx[-1] = xNodes[-1]-xNodes[-2]
        return dx
    
    
    def ComputeStretchedGridPoints(self, xCoords, xRefinement, location):      
        if location=='upstream':
            xNew = xCoords[0:-1].copy()
            dxMin = np.min(self.ComputeGridSpacing(xNew))
            dxRef = np.min(self.ComputeGridSpacing(xRefinement))
            while (dxMin-dxRef>0):
                newPoint =  xNew[-1] + (xRefinement[0]-xNew[-1])*0.5
                xNew = np.append(xNew, newPoint)
                dxMin = np.min(self.ComputeGridSpacing(xNew))
            xNew[-1] = xRefinement[0]
            # xNew = np.append(xNew, xRefinement[0])
        else:
            xNew = xCoords[1:].copy()
            dxMin = np.min(self.ComputeGridSpacing(xNew))
            dxRef = np.min(self.ComputeGridSpacing(xRefinement))
            while (dxMin-dxRef>0):
                newPoint =  xNew[0] - (+xNew[0] - xRefinement[-1])*0.5
                xNew = np.insert(xNew, 0, newPoint)
                dxMin = np.min(self.ComputeGridSpacing(xNew))
            xNew[0] = xRefinement[-1]
            # xNew = np.insert(xNew, 0, xRefinement[-1])
        return xNew
        
        
    
    def GenerateVirtualGeometry(self, xNodes):
        """
        Generate the virtual geometry consisting of halo nodes for boundary conditions

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.xNodes = xNodes
        self.nNodesHalo = self.nNodes+2
        self.xNodesVirt = np.zeros(self.nNodesHalo)
        self.xNodesVirt[1:-1] = self.xNodes
        self.xNodesVirt[0] = self.xNodes[0] - (xNodes[1]-xNodes[0])
        self.xNodesVirt[-1] = self.xNodes[-1] + (xNodes[-1]-xNodes[-2])
        self.areaReference = self.config.getAreaReference()
        self.dx = self.ComputeGridSpacing(self.xNodesVirt)
        
        if self.topology.lower()=='default':
            print("The simulation proceeds with default topology: constant area")
            self.areaTube = np.zeros_like(self.xNodesVirt)+self.areaReference
        elif self.topology.lower()=='nozzle':
            print(f"The simulation topology is: nozzle. Reading the coordinates from the nozzle file {self.config.getNozzleFilePath()}")
            self.areaTube = self.readNozzleFile(self.xNodesVirt, self.config.getNozzleFilePath())
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
            self.solution[name] = np.zeros(self.nNodesHalo)


    def InstantiateSolutionArraysConservatives(self):
        """
        Instantiate the containers for the solutions. The first dimension is space, the second is time.
        """
        self.solutionConsNames = ['u1', 'u2', 'u3']
        self.solutionCons = {}
        for name in self.solutionConsNames:
            self.solutionCons[name] = np.zeros(self.nNodesHalo)
        

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
            self.solution[name] = self.CopyInitialState(initialConditions[name][0], initialConditions[name][1])
    

    def InitialConditionsArrays(self, dictIn):
        """
        Initialize the conditions based on initial state, defined by arrays
        """
        dictIn['Energy'] = dictIn['Pressure'] / (self.gmma - 1) / dictIn['Density']
        for name in self.solutionNames:
            self.solution[name][:, 1:-1] = dictIn[name]
    
    
    def PlotGridGeometry(self, trueAspectRatio=True, pointsToJump=1, save_filename=None):
        """Plot the grid geometry. 1D tube, with thickness equal to the diameter of the tube
        """
        diameter = np.sqrt(4*self.areaTube/np.pi)
        yLower = np.zeros_like(self.xNodesVirt)-diameter/2
        yUpper = diameter/2
        
        plt.figure()
        plt.plot(self.xNodesVirt, yLower, 'k')
        plt.plot(self.xNodesVirt, yUpper, 'k')
        nPointsPic = 10
        for i in range(0, len(diameter), pointsToJump):
            plt.plot(np.zeros(nPointsPic)+self.xNodesVirt[i], np.linspace(yLower[i], yUpper[i], nPointsPic), '-k', lw=0.5)
        plt.xlabel(r'$x \ \rm{[m]}$')
        plt.ylabel(r'$r \ \rm{[m]}$')
        
        if trueAspectRatio:
            ax = plt.gca()
            ax.set_aspect('equal')
        
        if save_filename is not None:
            plt.savefig(save_filename + '.pdf', bbox_inches='tight')
        

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

    def SetBoundaryConditions(self):
        """
        Set the correct boundary condition type (`reflective`, `transparent`, or `periodic`)
        """
        if self.BCtype[0].lower()=='reflective':
            self.SetReflectiveBoundaryConditions('left')
        elif self.BCtype[0].lower()=='transparent':
            self.SetTransparentBoundaryConditions('left')
        elif self.BCtype[0].lower()=='periodic':
            self.SetPeriodicBoundaryConditions('left')
        elif self.BCtype[0].lower()=='inlet':
            self.SetInletBoundaryConditions('left')
        elif self.BCtype[0].lower()=='outlet':
            self.SetOutletBoundaryConditions('left')
        else:
            raise ValueError("Unknown boundary condition type on the left")
        
        if self.BCtype[1].lower()=='reflective':
            self.SetReflectiveBoundaryConditions('right')
        elif self.BCtype[1].lower()=='transparent':
            self.SetTransparentBoundaryConditions('right')
        elif self.BCtype[1].lower()=='periodic':
            self.SetPeriodicBoundaryConditions('right')
        elif self.BCtype[1].lower()=='outlet':
            self.SetOutletBoundaryConditions('right')
        elif self.BCtype[1].lower()=='inlet':
            self.SetInletBoundaryConditions('right')
        else:
            raise ValueError("Unknown boundary condition type on the right")
        
        # update also the conservative variable arrays based on what has been done on the primitive
        self.solutionCons['u1'], self.solutionCons['u2'], self.solutionCons['u3'] = (GetConservativesFromPrimitives(
            self.solution['Density'], self.solution['Velocity'], self.solution['Pressure'], self.fluid))


    def SetReflectiveBoundaryConditions(self, location):
        """
        Set reflective BC
        """
        if location=='left':
            self.solution['Density'][0] = self.solution['Density'][1]
            self.solution['Velocity'][0] = -self.solution['Velocity'][1]
            self.solution['Pressure'][0] = self.solution['Pressure'][1]
            self.solution['Energy'][0] = self.solution['Energy'][1]
        elif location=='right':
            self.solution['Density'][-1] = self.solution['Density'][-2]
            self.solution['Velocity'][-1] = -self.solution['Velocity'][-2]
            self.solution['Pressure'][-1] = self.solution['Pressure'][-2]
            self.solution['Energy'][-1] = self.solution['Energy'][-2]
        else:
            raise ValueError('Unknown location specified')
            
    

    def SetTransparentBoundaryConditions(self, location):
        """
        Set transparent BC
        """
        if location=='left':
            self.solution['Density'][0] = self.solution['Density'][1]
            self.solution['Velocity'][0] = self.solution['Velocity'][1]
            self.solution['Pressure'][0] = self.solution['Pressure'][1]
            self.solution['Energy'][0] = self.solution['Energy'][1]
        elif location=='right':
            self.solution['Density'][-1] = self.solution['Density'][-2]
            self.solution['Velocity'][-1] = self.solution['Velocity'][-2]
            self.solution['Pressure'][-1] = self.solution['Pressure'][-2]
            self.solution['Energy'][-1] = self.solution['Energy'][-2]
        else:
            raise ValueError('Unknown location specified')
        
        
        
    def SetPeriodicBoundaryConditions(self, location):
        """
        Set periodic BC
        """
        if location=='left':
            self.solution['Density'][0] = self.solution['Density'][-2]
            self.solution['Velocity'][0] = self.solution['Velocity'][-2]
            self.solution['Pressure'][0] = self.solution['Pressure'][-2]
            self.solution['Energy'][0] = self.solution['Energy'][-2]
        elif location=='right':
            self.solution['Density'][-1] = self.solution['Density'][1]
            self.solution['Velocity'][-1] = self.solution['Velocity'][1]
            self.solution['Pressure'][-1] = self.solution['Pressure'][1]
            self.solution['Energy'][-1] = self.solution['Energy'][1]
        else:
            raise ValueError('Unknown location specified')
    
    
    def SetInletBoundaryConditions(self, location):
        """
        Set periodic BC
        """
        # handle left and right extremities with the same code
        if location=='right':
            iHalo = -1
            iInternal = -2
        elif location=='left':
            iHalo = 0
            iInternal = 1
        else:
            raise ValueError('Unknown location specified')
        
        inletConditions = self.config.getInletConditions()
        totalPressure = inletConditions[0]
        totalTemperature = inletConditions[1]
        direction = inletConditions[2]
        
        # static pressure is the only info taken from the domain
        pressure = self.solution['Pressure'][iInternal]
        if isinstance(self.fluid, FluidIdeal):
            mach = self.fluid.ComputeMach_pt_p(totalPressure, pressure)
            temperature = self.fluid.ComputeTemperature_Tt_M(totalTemperature, mach)
            density = self.fluid.ComputeDensity_p_T(pressure, temperature)
        else:
            raise ValueError('At the moment inlet/outlet BC are implemented only for ideal gas')
            mach, temperature, density = self.fluid.SolveRealGasInlet(pressure, totalPressure, totalTemperature)
        
        soundSpeed = self.fluid.ComputeSoundSpeed_p_rho(pressure, density)
        velocity = mach*soundSpeed*direction
        energy = self.fluid.ComputeStaticEnergy_p_rho(pressure, density)
            
        self.solution['Density'][iHalo] = density
        self.solution['Velocity'][iHalo] = velocity
        self.solution['Pressure'][iHalo] = pressure
        self.solution['Energy'][iHalo] = energy
    
    
    def SetOutletBoundaryConditions(self, location):
        """
        Set periodic BC at time
        """
        # handle left and right extremities with the same code
        if location=='right':
            iHalo = -1
            iInternal = -2
        elif location=='left':
            iHalo = 0
            iInternal = 1
        else:
            raise ValueError('Unknown location specified')
            
        machOutlet = self.fluid.ComputeMach_u_p_rho(self.solution['Velocity'][iInternal], self.solution['Pressure'][iInternal], self.solution['Density'][iInternal])        
        if machOutlet<1:
            pressure = self.config.getOutletConditions() # the pressure is the information taken from outside
            velocity = self.solution['Velocity'][iInternal]
            density = self.solution['Density'][iInternal]
            energy = self.fluid.ComputeStaticEnergy_p_rho(pressure, density)        
            self.solution['Density'][iHalo] = density
            self.solution['Velocity'][iHalo] = velocity
            self.solution['Pressure'][iHalo] = pressure
            self.solution['Energy'][iHalo] = energy
        else:            
            self.SetTransparentBoundaryConditions(location) # the boundary is equivalent to a transparent condition
            
            


    def SolveSystem(self, high_order=False, limiter='Van Albada'):
        """
        Solve the equations explicitly in time (forward Euler) using a certain flux_method (`Godunov`, `Roe`, `WAF`). high_order
        specifies if applying or not high order reconstruction with limiters. At the moment only type one is working -> simply
        impose high_order=True
        """
        flux_method = self.config.getNumericalScheme()
        high_order = self.config.getMUSCLReconstruction()
        
        print()
        print("="*80)
        print(" "*33 + "START SOLVER")
        print("Numerical flux method: %s" %(flux_method))
        print("MUSCL reconstruction: %s" %high_order)
        print()

        # short aliases
        conservativeOld = self.solutionCons.copy()
        conservativeNew = self.solutionCons.copy()
        primitiveOld = self.solution.copy()
        primitiveNew = self.solution.copy()
        
        time = 0
        iTime = 1
        
        # time-steps loop
        while time < self.timeMax:
            dt = self.ComputeTimeStep(primitiveOld)
            if time + dt > self.timeMax:
                dt = self.timeMax - time
            newTime = time + dt
                    
            print(f"Iteration: {iTime}, Progress in Time {((newTime)/self.timeMax * 100):.3f} %")
            
            # compute fluxes on every internal interface
            flux = np.zeros((self.nNodes+1, 3))
            for iFace in range(flux.shape[0]):
                flux[iFace, :] = self.ComputeFluxVector(iFace, iFace+1, primitiveOld, dt, flux_method, high_order, limiter)
            
            if self.topology.lower()=='nozzle':
                source = self.ComputeSourceTerms(primitiveOld)
            else:
                source = np.zeros((self.nNodesHalo,3))
            
            # Update the conservatives for every physical cell element (vectorized version)
            conservativeNew['u1'][1:-1] = conservativeOld['u1'][1:-1] + dt/self.dx[1:-1] * ((flux[0:-1, 0] - flux[1:, 0]) + source[1:-1, 0]*self.dx[1:-1])
            conservativeNew['u2'][1:-1] = conservativeOld['u2'][1:-1] + dt/self.dx[1:-1] * ((flux[0:-1, 1] - flux[1:, 1]) + source[1:-1, 1]*self.dx[1:-1])
            conservativeNew['u3'][1:-1] = conservativeOld['u3'][1:-1] + dt/self.dx[1:-1] * ((flux[0:-1, 2] - flux[1:, 2]) + source[1:-1, 2]*self.dx[1:-1])

            # update the primitives
            self.solution['Density'][1:-1], self.solution['Velocity'][1:-1], self.solution['Pressure'][1:-1], self.solution['Energy'][1:-1] = \
                GetPrimitivesFromConservatives(conservativeNew['u1'][1:-1], conservativeNew['u2'][1:-1], conservativeNew['u3'][1:-1], self.fluid)
            
            self.checkSimulationStatus()
            
            # set boundary conditions to update the ghost points for the new iteration
            self.SetBoundaryConditions()
            
            self.WriteSolution(iTime, newTime)
            
            time += dt  
            iTime += 1
            
        print(" "*34 + "END SOLVER")
        print("="*80)
    
    
    def ComputeTimeStep(self, primitive):
        velocity = primitive['Velocity'][1:-1]
        speedOfSound = self.fluid.ComputeSoundSpeed_p_rho(primitive['Pressure'][1:-1], primitive['Density'][1:-1])
        dtMax = np.min(self.dx[1:-1] * self.cflMax / (np.abs(velocity)+speedOfSound))
        return dtMax
    
    
    def WriteSolution(self, it, time):        
        full_path = self.subFolder + '/' + 'step_%06i.pik' %(it)
        outputResults = {'Time': time, 'Iteration Counter': it, 'X Coords': self.xNodesVirt, 'Primitive': self.solution, 'Fluid': self.fluid}
        with open(full_path, 'wb') as file:
            pickle.dump(outputResults, file)
        
    
    
    def removeUnusedSpace(self, itMax):
        for key in self.solutionCons.keys():
            self.solutionCons[key] = self.solutionCons[key][:, :itMax+1]
        
        for key in self.solution.keys():
            self.solution[key] = self.solution[key][:, :itMax+1]

        self.timeVec = self.timeVec[:itMax+1]
        self.nTime = len(self.timeVec)
    
    
    def ComputeSourceTerms(self, primitive):
        """Compute source terms related to area variations along the tube due to a nozzle. Source terms taken from 'On the numerical simulation
        of non-classical quasi-1D steady nozzle flows: Capturing sonic shocks' by Vimercati and Guardone.

        Args:
            it (int): time step index

        Returns:
            np.ndarray: source terms arrays (nPoints, 3)
        """
        totalEnergy = primitive['Energy'][:] + 0.5*primitive['Velocity']**2
        source_terms = np.zeros((self.nNodesHalo,3))
        source_terms[:,0] = - primitive['Density'] * primitive['Velocity']*self.dAreaTude_dx/self.areaTube
        source_terms[:,1] = - (primitive['Density'] * primitive['Velocity']**2)*self.dAreaTude_dx/self.areaTube
        source_terms[:,2] = - primitive['Velocity'] *(primitive['Density']*totalEnergy + primitive['Pressure'])*self.dAreaTude_dx/self.areaTube
        return source_terms
    
    
    def checkSimulationStatus(self):
        """
        Check if nans or infs are detected and in that case stop the simulation and provide explanation
        """
        if np.any(np.isnan(self.solution['Density'])) or np.any(np.isinf(self.solution['Density'])) or \
            np.any(np.isnan(self.solution['Pressure'])) or np.any(np.isinf(self.solution['Pressure'])):
            print()
            print()
            print("######################  SIMULATION DIVERGED ############################")
            print('NaNs detected in density. Simulation stopped.')
            cfl = self.ComputeMaxCFL() # use the previous time step to compute where the solution had CFL related problems
            dum = it-1
            while np.any(np.isnan(cfl)):
                dum -=1
                cfl = self.ComputeMaxCFL(dum)
            print("Maximum CFL number found: %.3f" %(np.max(cfl)))
            print("At location x: %.3f [m], at time: %.3e [s]" %(self.xNodesVirt[np.argmax(cfl)], self.timeVec[dum]))
            print("Visualize the plot to understand critical locations, and decrease CFL_MAX input setting.")
            print("###############################  EXIT ##################################")
            print()
            
            plt.figure()
            plt.plot(self.xNodes, cfl)
            plt.xlabel('x [m]')
            plt.ylabel('CFL [-]')
            plt.grid(alpha=.3)
            plt.show()
            sys.exit()
    
    
    def ComputeMaxCFL(self, it):
        pressure = self.solution['Pressure'][1:-1,it]
        density = self.solution['Density'][1:-1,it]
        velocity = self.solution['Velocity'][1:-1,it]
        dx = self.dx[1:-1]
        soundSpeed = np.zeros_like(pressure)
        for i in range(len(soundSpeed)):
            soundSpeed = self.fluid.ComputeSoundSpeed_p_rho(pressure[i], density[i])
        cfl = (np.abs(velocity)+soundSpeed)*self.dt/dx
        return cfl
        


    def ComputeFluxVector(self, il, ir, primitive, dt, flux_method, high_order, limiter):
        """
        Compute the flux vector at the interface between grid points `il` and `ir`, using a certain `flux_method`.
        """
        
        # flow reconstruction if high_order=True
        if (high_order and il>2 and ir<self.nNodesHalo-2):
            rhoL, uL, pL, rhoR, uR, pR = self.MUSCL_VanAlbada_Reconstruction(il, ir, it)  # currently working
            # rhoL, uL, pL, rhoR, uR, pR = self.MUSCL_Reconstruction(il, ir, it, limiter)  # at the moment creates oscillations with every limiter for some reasons
        else:
            rhoL = primitive['Density'][il]
            rhoR = primitive['Density'][ir]
            uL = primitive['Velocity'][il]
            uR = primitive['Velocity'][ir]
            pL = primitive['Pressure'][il]
            pR = primitive['Pressure'][ir]            
        
        # flux calculation
        if flux_method.lower()=='godunov':
            if self.fluid_model!='ideal':
                raise ValueError('Godunov scheme is available only for ideal gas model')
            nx, nt = 51, 51
            x = np.linspace(-self.dx/2, self.dx/2, nx)
            t = np.linspace(0, dt, nt)
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
            t = np.linspace(0, dt, nt)
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
                rhoLL, rhoRR = primitive['Density'][il-1], primitive['Density'][ir+1]
                uLL, uRR = primitive['Velocity'][il-1], primitive['Velocity'][ir+1]
                pLL, pRR = primitive['Pressure'][il-1], primitive['Pressure'][ir+1]
            else: # no extrapolation for extreme cells over the halo nodes
                rhoLL, rhoRR = rhoL, rhoR
                uLL, uRR = uL, uR
                pLL, pRR = pL, pR
            mhck = MusclHancock(rhoLL, rhoL, rhoR, rhoRR, uLL, uL, uR, uRR, pLL, pL, pR, pRR, self.dx)
            mhck.ReconstructInterfaceValues()
            mhck.EvolveInterfaceValues(self.dx, dt)
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


    # def ShowAnimation(self):
    #     """
    #     Show animation of the results at all time instants
    #     """
    #     ni, nt = self.solution['Density'].shape
    #     mach = np.zeros((ni, nt))
    #     for i in range(ni):
    #         for t in range(nt):
    #             mach[i,t] = self.solution['Velocity'][i,t]/self.fluid.ComputeSoundSpeed_p_rho(self.solution['Pressure'][i,t], self.solution['Density'][i,t])
    #     def plot_limits(f, extension=0.05):
    #         max = f.max()
    #         min = f.min()
    #         left = min-(max-min)*extension
    #         right = max+(max-min)*extension
    #         return left, right
        
    #     if self.config.showAnimation():
    #         fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    #         density_limits = plot_limits(self.solution['Density'])
    #         velocity_limits = plot_limits(self.solution['Velocity'])
    #         pressure_limits = plot_limits(self.solution['Pressure'])
    #         energy_limits = plot_limits(self.solution['Energy'])
    #         mach_limitis = plot_limits(mach)
    #         for it in range(self.nTime):
    #             for row in ax:
    #                 for col in row:
    #                     col.cla()
    #             ax[0, 0].plot(self.xNodesVirt, self.solution['Density'][:, it], '-C0o', ms=2)
    #             ax[0, 0].set_ylabel(r'Density [kg/m3]')
    #             ax[0, 0].set_ylim(density_limits)

    #             ax[0, 1].plot(self.xNodesVirt, self.solution['Velocity'][:, it], '-C1o', ms=2)
    #             ax[0, 1].set_ylabel(r'Velocity [m/s]')
    #             ax[0, 1].set_ylim(velocity_limits)

    #             ax[1, 0].plot(self.xNodesVirt, self.solution['Pressure'][:, it], '-C2o', ms=2)
    #             ax[1, 0].set_ylabel(r'Pressure [Pa]')
    #             ax[1, 0].set_ylim(pressure_limits)

    #             # ax[1, 1].plot(self.xNodesVirt, self.solution['Energy'][:, it], '-C3o', ms=2)
    #             # ax[1, 1].set_ylabel(r'Energy')
    #             # ax[1, 1].set_ylim(energy_limits)
    #             ax[1, 1].plot(self.xNodesVirt, mach[:, it], '-C3o', ms=2)
    #             ax[1, 1].set_ylabel(r'Mach [-]')
    #             ax[1, 1].set_ylim(mach_limitis)

    #             fig.suptitle('Time %.3e [s]' % self.timeVec[it])

    #             for row in ax:
    #                 for col in row:
    #                     col.set_xlabel('x')
    #                     col.grid(alpha=.3)
    #             plt.pause(1e-3)
    #         plt.show()
    

    # def ShowContourAnimation(self):
    #     """
    #     Show contour animation of the results for all time instants
    #     """
    #     xgrid = np.reshape(self.xNodes, (len(self.xNodes), 1))
    #     xgrid = np.hstack((xgrid, xgrid))
    #     ygrid = np.reshape(np.zeros_like(self.xNodes), (len(self.xNodes), 1))
    #     length = self.xNodes[-1]-self.xNodes[0]
    #     AR = 4
    #     height = length/AR
    #     ygrid = np.hstack((ygrid, ygrid+height))

    #     def convert2D(f):
    #         fgrid = np.reshape(f, (len(self.xNodes), 1))
    #         fgrid = np.hstack((fgrid, fgrid))
    #         return fgrid

    #     fig, ax = plt.subplots(2, 2, figsize=(8, 5))
    #     for it in range(self.nTime):
    #         for row in ax:
    #             for col in row:
    #                 col.cla()
    #         rhoGrid = convert2D(self.solution['Density'][1:-1, it])
    #         uGrid = convert2D(self.solution['Velocity'][1:-1, it])
    #         pGrid = convert2D(self.solution['Pressure'][1:-1, it])
    #         eGrid = convert2D(self.solution['Energy'][1:-1, it])
           
    #         ax[0, 0].contourf(xgrid, ygrid, rhoGrid, cmap='jet', levels=20)
    #         ax[0, 0].set_title(r'Density')

    #         ax[0, 1].contourf(xgrid, ygrid, uGrid, cmap='jet', levels=20)
    #         ax[0, 1].set_title(r'Velocity')

    #         ax[1, 0].contourf(xgrid, ygrid, pGrid, cmap='jet', levels=20)
    #         ax[1, 0].set_title(r'Pressure')

    #         ax[1, 1].contourf(xgrid, ygrid, eGrid, cmap='jet', levels=20)
    #         ax[1, 1].set_title(r'Energy')

    #         fig.suptitle('Time %.3f' % self.timeVec[it])

    #         for row in ax:
    #             for col in row:
    #                 col.set_aspect('equal')
    #         plt.pause(1e-3)
    

    def SaveSolution(self):
        """
        Save the full object as a pickle for later use
        """
        folder_name = self.config.getOutputFolder()
        os.makedirs(folder_name, exist_ok=True)
        file_name = self.config.getOutputFileName()
        full_path = folder_name+'/'+file_name+'_NX_%i_TMAX_%.6f.pik' %(self.nNodes, self.timeMax)
        with open(full_path, 'wb') as file:
            pickle.dump(self, file)
        print('Pickle object with full solution saved to ' + full_path + ' !')


    def SaveNodeSolutionToCSV(self, iNode, timeInstants, folder_name, file_name):
        """
        Save the array of fluid flow quantities (P,T,s,Mach,Gamma) from the solution to a CSV file.
        """
        file_path = folder_name + '/' + file_name + '.dat'

        pressure_data = self.solution['Pressure'][iNode, :]  # Extract the pressure data (1D array)
        density_data = self.solution['Density'][iNode, :]  # Extract the density data (1D array)
        temperature_data = self.fluid.ComputeTemperature_p_rho(pressure_data, density_data)
        entropy_data = self.fluid.ComputeEntropy_p_rho(pressure_data,density_data)
        fundDerGasDynamics_data = self.fluid.ComputeFunDerGamma_p_rho(pressure_data,density_data)
        compressibilityFactor_data = self.fluid.ComputeComprFactorZ_p_rho(pressure_data,density_data)

        with open(file_path, 'w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            for value in range(len(timeInstants)):
                writer.writerow([timeInstants[value], pressure_data[value], temperature_data[value], density_data[value],
                                 entropy_data[value], fundDerGasDynamics_data[value], compressibilityFactor_data[value]])

        print(f"Fluid flow quantities (P,T,D,s,Gamma,Z) saved to {file_path}!")

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
        
        # Linear interpolation with external filling set to area Reference (=Tube area)
        areaReference = self.config.getAreaReference()
        interpolatedNozzleArea = np.interp(xTube, nozzleX, nozzleArea, left=areaReference, right=areaReference)
    
        print(f"The reference tube area is: {areaReference:.6f} [m2].")
        print(f"The nozzle throat area is {interpolatedNozzleArea.min():.6f} [m2].")
        print(f"The nozzle maximum area is {interpolatedNozzleArea.max():.6f} [m2].")
        print(f"The area ratio between nozzle throat and exit section is {interpolatedNozzleArea.min()/interpolatedNozzleArea[-1]:.6f}.")
        print(f"The area ratio between nozzle throat and tube is {interpolatedNozzleArea.min()/areaReference:.6f}.")
        print(f"If this is not correct, modify the REFERENCE_AREA setting in the geometry section of the input file to the correct value for the tube area, or modify the nozzle csv file to be consistent with the tube area.")
        
        return interpolatedNozzleArea
        























