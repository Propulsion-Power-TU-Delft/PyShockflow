import matplotlib.pyplot as plt
import numpy as np
from PyShockTube.shock_tube import ShockTube

"""
INPUT PARAMETERS FOR THE SHOCK-TUBE PROBLEM
"""
LENGTH = 1
NX = 100
TIME_MAX = 0.25  # to simulate two reflections
RHOL, RHOR = 1.0, 0.125
UL, UR = 0.0, 0.0
PL, PR = 1.0, 0.1
NUMERICAL_SCHEME = 'roe'
BOUNDARY_CONDITIONS = 'reflective'
FLUID = 'air'
FLUID_MODEL = 'ideal'
FLUID_GAMMA = 1.4
CFL_MAX = [0.01, 0.05, 0.10, 0.2, 0.5]
HIGH_ORDER = True




""" 
Solution Driver
"""
x = np.linspace(0, LENGTH, NX)
dx = x[1]-x[0]
Ssound = np.sqrt(np.array([FLUID_GAMMA*PL/RHOL, 
                           FLUID_GAMMA*PR/RHOR]))
Smax = np.max(Ssound)

for i in range(len(CFL_MAX)):
    dtMax = CFL_MAX[i]* dx / Smax
    nt = int(TIME_MAX/dtMax)
    t = np.linspace(0, TIME_MAX, nt)

    tube = ShockTube(x, t, FLUID, FLUID_MODEL, FLUID_GAMMA)
    inCondDict = {'Density': np.array([RHOL, RHOR]), 'Velocity': np.array([UL, UR]), 'Pressure': np.array([PL, PR])}
    tube.InstantiateSolutionArrays()
    tube.InstantiateSolutionArraysConservatives()
    tube.InitialConditionsLeftRight(inCondDict)
    tube.SetBoundaryConditions(BOUNDARY_CONDITIONS, 0)
    tube.SolveSystem(flux_method=NUMERICAL_SCHEME, high_order=HIGH_ORDER)
    tube.SaveSolution(folder_name='solutions', file_name='SodsTest_%.2f' %CFL_MAX[i])

plt.show()