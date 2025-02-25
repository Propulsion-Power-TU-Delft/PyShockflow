import matplotlib.pyplot as plt
import numpy as np
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
