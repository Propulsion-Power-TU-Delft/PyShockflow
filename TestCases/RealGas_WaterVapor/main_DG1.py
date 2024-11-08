import matplotlib.pyplot as plt
import numpy as np
from PyShockTube.shock_tube import ShockTube
import CoolProp.CoolProp as CP

"""
INPUT PARAMETERS FOR THE SHOCK-TUBE PROBLEM
"""
LENGTH = 1
NX = 250
UL, UR = 0.0, 0.0
Rhocr = 322
RHOL = 1.818*Rhocr
RHOR = 0.275*Rhocr
Pcr = 22.064e6
PL, PR = 3*Pcr, 0.575*Pcr
TIME_MAX = 0.15*LENGTH*(Rhocr/Pcr)**(0.5)
CFLmax = 0.95  
FLUID = 'water'
FLUID_MODEL = 'real'
NUMERICAL_SCHEME = 'roe'
BOUNDARY_CONDITIONS = 'reflective'


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
tube.SetBoundaryConditions(BOUNDARY_CONDITIONS, 0)
tube.SolveSystem(flux_method=NUMERICAL_SCHEME)
tube.SaveSolution(folder_name='solutions', file_name='DG1')
# tube.ShowAnimation()




plt.show()