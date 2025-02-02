import matplotlib.pyplot as plt
import numpy as np
from PyShockTube.shock_tube import ShockTube
import CoolProp.CoolProp as CP

"""
INPUT PARAMETERS FOR THE SHOCK-TUBE PROBLEM
"""
LENGTH = 3.6                            # minimal length to include the four pressure probes (PT1/2/3/4)
NX = 1000                                # number of discretization points
TIME_MAX = 0.06                         # adjust as needed
# Gamma = -0.219
TL, TR = 640.65, 573.15                 # left and right initial temperatures
UL, UR = 0.0, 0.0                       # left and right initial velocities
PL, PR = 8.9E+05, 0.02E+05              # left and right initial pressures
# Gamma = -0.016
# TL, TR = 645.65, 573.15                 # left and right initial temperatures
# UL, UR = 0.0, 0.0                       # left and right initial velocities
# PL, PR = 9.3E+05, 0.02E+05              # left and right initial pressures
NUMERICAL_SCHEME = 'roe'                # the only scheme working with real gas
BOUNDARY_CONDITIONS = 'transparent'     # to mimic the effect of the discharge nozzle
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
tube.SolveSystem(flux_method=NUMERICAL_SCHEME, high_order=False)
tube.SaveSolution(folder_name='solutions', file_name='ASTER_D6_%.2f' %TIME_MAX)
# tube.ShowAnimation()

# Plot the trend of pressure over time as measured by the Kulite sensors
# Positions of the Kulite sensors of the ASTER from the diaphragm
iNode_Kulite1 = int( NX* ((LENGTH/2) - 1.2) / LENGTH ) # Kulite 1 (1200 mm from diaphragm)
iNode_Kulite2 = int( NX* ((LENGTH/2) - 1.6) / LENGTH ) # Kulite 2 (1600 mm from diaphragm)
iNode_Kulite3 = int( NX* ((LENGTH/2) - 1.8) / LENGTH ) # Kulite 3 (1800 mm from diaphragm)
# print(iNode_Kulite1, iNode_Kulite2, iNode_Kulite3)

tube.PlotNodeSolution(iNode_Kulite1, folder_name = 'solutions', file_name = 'solution_at_Kulite1')
tube.PlotNodeSolution(iNode_Kulite2, folder_name = 'solutions', file_name = 'solution_at_Kulite2')
tube.PlotNodeSolution(iNode_Kulite3, folder_name = 'solutions', file_name = 'solution_at_Kulite3')
tube.SaveNodeSolutionToCSV(iNode_Kulite1, t, folder_name='solutions', file_name='pressure_trend_at_Kulite1')
tube.SaveNodeSolutionToCSV(iNode_Kulite2, t, folder_name='solutions', file_name='pressure_trend_at_Kulite2')
tube.SaveNodeSolutionToCSV(iNode_Kulite3, t, folder_name='solutions', file_name='pressure_trend_at_Kulite3')

# plt.show()
