from PyShockflow.shock_tube import ShockTube
from PyShockflow.config import Config

config = Config('input_1stOrder.ini')
tube = ShockTube(config)
tube.InstantiateSolutionArrays()
tube.InstantiateSolutionArraysConservatives()
tube.InitialConditionsLeftRight()
tube.SetBoundaryConditions()
tube.SolveSystem()
