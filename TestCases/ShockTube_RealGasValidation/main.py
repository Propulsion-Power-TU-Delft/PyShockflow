import matplotlib.pyplot as plt
import numpy as np
from NumericalCodes.shock_tube import ShockTube
import CoolProp.CoolProp as CP

"""
INPUT PARAMETERS FOR THE SHOCK-TUBE PROBLEM
"""
LENGTH = 10
NX = 500
RHOL, RHOR = 348.8, 3.488
UL, UR = 0.0, 0.0
TL, TR = 892.67, 1116.89
TIME_MAX = 0.00275
# TIME_MAX = 0.000275*3
CFLmax = 0.85  # conservative CFL
FLUID = 'co2'
FLUID_MODEL = 'real'


""" 
Solution Driver
"""
PL = CP.PropsSI('P', 'D', RHOL, 'T', TL, FLUID)
PR = CP.PropsSI('P', 'D', RHOR, 'T', TR, FLUID)
PL, PR = 73760000, 737600
x = np.linspace(0, LENGTH, NX)
dx = x[1]-x[0]
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
tube.SaveSolution(folder_name='solutions', file_name='co2_entropyfix_%i' %NX)
tube.ShowAnimation()




plt.show()