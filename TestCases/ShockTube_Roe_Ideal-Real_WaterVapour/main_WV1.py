import matplotlib.pyplot as plt
import numpy as np
from NumericalCodes.shock_tube import ShockTube
import CoolProp.CoolProp as CP

"""
INPUT PARAMETERS FOR THE SHOCK-TUBE PROBLEM
"""
LENGTH = 1
NX = 350
UL, UR = 0.0, 0.0
Rhocr = 322
RHOL = 1.01*Rhocr
RHOR = 0.59400*Rhocr
Pcr = 22.064e6
PL, PR = 1.6077*Pcr, 0.8957*Pcr
TIME_MAX = 0.2*LENGTH*(Rhocr/Pcr)**(0.5)
CFLmax = 0.5  # conservative CFL
FLUID = 'water'
FLUID_MODEL = 'real'


""" 
Solution Driver
"""
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
tube.SaveSolution(folder_name='solutions', file_name='WV1')
# tube.ShowAnimation()




plt.show()