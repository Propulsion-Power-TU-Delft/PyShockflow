from PyShockTube.shock_tube import ShockTube
from PyShockTube.config import Config

pressureList = [94]
configList = ['input_%ikPa.ini' %(pressure) for pressure in pressureList]

for configFile in configList:
    config = Config(configFile)
    tube = ShockTube(config)
    tube.InstantiateSolutionArrays()
    tube.InstantiateSolutionArraysConservatives()
    tube.InitialConditionsLeftRight()
    tube.SetBoundaryConditions(0)
    tube.SolveSystem()
    tube.SaveSolution()
    tube.ShowAnimation()
