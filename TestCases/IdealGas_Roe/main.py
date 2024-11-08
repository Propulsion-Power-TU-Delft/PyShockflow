import matplotlib.pyplot as plt
import numpy as np
from PyShockTube.shock_tube import ShockTube

"""
INPUT PARAMETERS FOR THE SHOCK-TUBE PROBLEM
"""
LENGTH = 1
NX = 100
TIME_MAX = 2  # to simulate two reflections
RHOL, RHOR = 1.0, 0.125
UL, UR = 0.0, 0.0
PL, PR = 1.0, 0.1
NUMERICAL_SCHEME = 'roe'
BOUNDARY_CONDITIONS = 'reflective'
FLUID = 'air'
FLUID_MODEL = 'ideal'
FLUID_GAMMA = 1.4




""" 
Solution Driver
"""
x = np.linspace(0, LENGTH, NX)
dx = x[1]-x[0]
Smax = np.sqrt(1.4*np.max([PL, PR])/np.min([RHOL, RHOR]))+np.max([UL, UR])  # brutal approximation max eigenvalue
CFLmax = 0.9  # conservative CFL
dtMax = CFLmax* dx / Smax
nt = int(TIME_MAX/dtMax)
t = np.linspace(0, TIME_MAX, nt)

tube = ShockTube(x, t, FLUID, FLUID_MODEL, FLUID_GAMMA)
inCondDict = {'Density': np.array([RHOL, RHOR]), 'Velocity': np.array([UL, UR]), 'Pressure': np.array([PL, PR])}
tube.InstantiateSolutionArrays()
tube.InstantiateSolutionArraysConservatives()
tube.InitialConditionsLeftRight(inCondDict)
tube.SetBoundaryConditions(BOUNDARY_CONDITIONS, 0)
tube.SolveSystem(flux_method=NUMERICAL_SCHEME)
tube.SaveSolution(folder_name='solutions', file_name='SodsTest_tMax_%.2f' %TIME_MAX)
tube.ShowAnimation()




plt.show()