import matplotlib.pyplot as plt
import numpy as np
from NumericalCodes.shock_tube import ShockTube
import CoolProp.CoolProp as CP

"""
INPUT PARAMETERS FOR THE SHOCK-TUBE PROBLEM
"""
LENGTH = 5
NX = 150
TIME_MAX = 0.05  # to simulate two reflections
UL, UR = 0.0, 0.0
TL, TR = 500, 500
PL, PR = 50e5, 10e5
CFLmax = 0.9  # conservative CFL
FLUID = 'R134a'
FLUID_MODEL = 'real'



""" 
Solution Driver
"""
x = np.linspace(0, LENGTH, NX)
dx = x[1]-x[0]
RHOL = CP.PropsSI('D', 'P', PL, 'T', TL, FLUID)
RHOR = CP.PropsSI('D', 'P', PR, 'T', TR, FLUID)
Smax = np.sqrt(1.4*np.max([PL, PR])/np.min([RHOL, RHOR]))+np.max([UL, UR])  # brutal approximation max eigenvalue
dtMax = CFLmax* dx / Smax
nt = int(TIME_MAX/dtMax)
t = np.linspace(0, TIME_MAX, nt)

tube = ShockTube(x, t, FLUID, FLUID_MODEL)
inCondDict = {'Density': np.array([RHOL, RHOR]), 'Velocity': np.array([UL, UR]), 'Pressure': np.array([PL, PR])}
tube.InstantiateSolutionArrays()
tube.InstantiateSolutionArraysConservatives()
tube.InitialConditionsLeftRight(inCondDict)
tube.SetBoundaryConditions('reflective', 0)

tube.SolveSystem(flux_method='Roe')
tube.SaveSolution(folder_name='solutions', file_name='real_gas_%.3f' %TIME_MAX)
tube.ShowAnimation()




plt.show()