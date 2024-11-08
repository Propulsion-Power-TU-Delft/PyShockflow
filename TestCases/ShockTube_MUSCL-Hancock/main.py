import matplotlib.pyplot as plt
import numpy as np
from NumericalCodes.shock_tube import ShockTube
import time

"""
INPUT PARAMETERS FOR THE SHOCK-TUBE PROBLEM
"""
LENGTH = 1
NX = 100
TIME_MAX = 0.25  # to simulate two reflections
RHOL, RHOR = 1.0, 0.125
UL, UR = 0.0, 0.0
PL, PR = 1.0, 0.1




""" 
Solution Driver
"""
x = np.linspace(0, LENGTH, NX)
dx = x[1]-x[0]
Smax = np.sqrt(1.4*np.max([PL, PR])/np.min([RHOL, RHOR]))+np.max([UL, UR])  # brutal approximation max eigenvalue
CFLmax = 0.1  # conservative CFL
dtMax = CFLmax* dx / Smax
nt = int(TIME_MAX/dtMax)
t = np.linspace(0, TIME_MAX, nt)

tube = ShockTube(x, t)
inCondDict = {'Density': np.array([RHOL, RHOR]), 'Velocity': np.array([UL, UR]), 'Pressure': np.array([PL, PR])}
tube.InstantiateSolutionArrays()
tube.InstantiateSolutionArraysConservatives()
tube.InitialConditionsLeftRight(inCondDict)
tube.SetBoundaryConditions('reflective', 0)

start = time.time()
tube.SolveSystem(flux_method='MUSCL-Hancock')
end = time.time()
print('Time: %.2f' %(end-start))
tube.SaveSolution(folder_name='solutions', file_name='SodsTest_tMax_%.2f' %TIME_MAX)
tube.ShowAnimation()


plt.show()