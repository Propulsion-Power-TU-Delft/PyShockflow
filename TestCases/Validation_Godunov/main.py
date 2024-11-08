import matplotlib.pyplot as plt
import numpy as np
from PyShockTube.shock_tube import ShockTube

"""
INPUT PARAMETERS FOR THE SHOCK-TUBE PROBLEM
"""
LENGTH = 1
NX = 250
FLUID = 'air'
FLUID_MODEL = 'ideal'
FLUID_GAMMA = 1.4

# Page 129 of Riemann Solvers and numerical methods for fluid dynamics by Toro et al.
initialCond = {'Test1': [1.0, 0.125, 0.0, 0.0, 1.0, 0.1],
               'Test2': [1.0, 1.0, -2.0, 2.0, 0.4, 0.4],
               'Test3': [1.0, 1.0, 0.0, 0.0, 1000.0, 0.01],
               'Test4': [1.0, 1.0, 0.0, 0.0, 0.01, 100.0],
               'Test5': [5.99924, 5.99242, 19.5975, -6.19633, 460.894, 46.0950]}

timeVectors = {'Test1': 0.25,
               'Test2': 0.15,
               'Test3': 0.012,
               'Test4': 0.035,
               'Test5': 0.035}

for key in initialCond.keys():
    inCond = initialCond[key]
    TIME_MAX = timeVectors[key]
    x = np.linspace(0, LENGTH, NX)
    dx = x[1]-x[0]
    PL, PR = inCond[4], inCond[5]
    RHOL, RHOR = inCond[0], inCond[1]
    UL, UR = inCond[2], inCond[3]
    approxSpeed = np.sqrt(1.4*np.max([PL, PR])/np.min([RHOL, RHOR]))  # brutal approximation max eigenvalue
    CFLmax = 0.1  # conservative CFL
    dtMax = CFLmax* dx / approxSpeed
    nt = int(TIME_MAX/dtMax)
    t = np.linspace(0, TIME_MAX, nt)

    #driver
    tube = ShockTube(x, t, FLUID, FLUID_MODEL, FLUID_GAMMA)
    inCondDict = {'Density': np.array([RHOL, RHOR]), 'Velocity': np.array([UL, UR]), 'Pressure': np.array([PL, PR])}
    tube.InstantiateSolutionArrays()
    tube.InstantiateSolutionArraysConservatives()
    tube.InitialConditionsLeftRight(inCondDict)
    tube.SetBoundaryConditions('reflective', 0)

    tube.SolveSystem(flux_method='Godunov')
    tube.SaveSolution(folder_name='solutions', file_name=key)

    # tube.ShowAnimation()




plt.show()