import matplotlib.pyplot as plt
import numpy as np
from PyShockTube.shock_tube import ShockTube
import CoolProp.CoolProp as CP

"""
INPUT PARAMETERS FOR THE SHOCK-TUBE PROBLEM
"""
LENGTH = 3.6                            # minimal length to include the four pressure probes (PT1/2/3/4)
NX = 250                                # number of discretization points
TIME_MAX = 0.01                         # adjust as needed
TL, TR = 273.15+349.8, 273.15+349.8     # left and right initial temperatures
UL, UR = 0.0, 0.0                       # left and right initial velocities
PL, PR = 2.66e5, 1e4                    # left and right initial pressures
NUMERICAL_SCHEME = 'roe'                # the only scheme working with real gas
BOUNDARY_CONDITIONS = 'transparent'     # to mimick the effect of the discharge nozzle
FLUID = 'D6'                            # fluid name following coolprop
FLUID_MODEL = 'real'                    # real gas model




""" 
Solution Driver
"""
RHOR = CP.PropsSI('D', 'T', TR, 'P', PR, FLUID)
RHOL = CP.PropsSI('D', 'T', TL, 'P', PL, FLUID)
x = np.linspace(0, LENGTH, NX)
dx = x[1]-x[0]
Smax = np.sqrt(1.4*np.max([PL, PR])/np.min([RHOL, RHOR]))+np.max([UL, UR])  # brutal approximation max eigenvalue, modify if needed
CFLmax = 0.9  # reduce if needed
dtMax = CFLmax* dx / Smax
nt = int(TIME_MAX/dtMax)
t = np.linspace(0, TIME_MAX, nt)
tube = ShockTube(x, t, FLUID, FLUID_MODEL)
inCondDict = {'Density': np.array([RHOL, RHOR]), 'Velocity': np.array([UL, UR]), 'Pressure': np.array([PL, PR])}
tube.InstantiateSolutionArrays()
tube.InstantiateSolutionArraysConservatives()
tube.InitialConditionsLeftRight(inCondDict)
tube.SetBoundaryConditions(BOUNDARY_CONDITIONS, 0)
tube.SolveSystem(flux_method=NUMERICAL_SCHEME)
tube.SaveSolution(folder_name='solutions', file_name='SodsTest_tMax_%.2f' %TIME_MAX)
tube.ShowAnimation()




plt.show()