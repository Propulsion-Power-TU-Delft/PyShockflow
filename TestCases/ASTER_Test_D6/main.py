from PyShockTube.shock_tube import ShockTube
from PyShockTube.config import Config

config = Config('input.ini')
tube = ShockTube(config)
tube.InstantiateSolutionArrays()
tube.InstantiateSolutionArraysConservatives()
tube.InitialConditionsLeftRight()
tube.SetBoundaryConditions(0)
tube.SolveSystem()
tube.SaveSolution()
tube.ShowAnimation()




# COPIED FROM main_deprecated.py

# Plot the trend of pressure over time as measured by the Kulite sensors
# Positions of the Kulite sensors of the ASTER from the diaphragm
NX = config.getNumberOfPoints()
LENGTH = config.getLength()
t = tube.timeVec

iNode_Kulite1 = int( NX* ((LENGTH/2) - 1.2) / LENGTH ) # Kulite 1 (1200 mm from diaphragm)
iNode_Kulite2 = int( NX* ((LENGTH/2) - 1.4) / LENGTH ) # Kulite 2 (1600 mm from diaphragm)
iNode_Kulite3 = int( NX* ((LENGTH/2) - 1.8) / LENGTH ) # Kulite 3 (1800 mm from diaphragm)
# print(iNode_Kulite1, iNode_Kulite2, iNode_Kulite3)

tube.PlotNodeSolution(iNode_Kulite1, folder_name = 'solutions', file_name = 'solution_at_Kulite1')
tube.PlotNodeSolution(iNode_Kulite2, folder_name = 'solutions', file_name = 'solution_at_Kulite2')
tube.PlotNodeSolution(iNode_Kulite3, folder_name = 'solutions', file_name = 'solution_at_Kulite3')
tube.SaveNodeSolutionToCSV(iNode_Kulite1, t, folder_name='solutions', file_name='PT_trend_at_Kulite1')
tube.SaveNodeSolutionToCSV(iNode_Kulite2, t, folder_name='solutions', file_name='PT_trend_at_Kulite2')
tube.SaveNodeSolutionToCSV(iNode_Kulite3, t, folder_name='solutions', file_name='PT_trend_at_Kulite3')