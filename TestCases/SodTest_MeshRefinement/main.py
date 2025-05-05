from PyShockflow.shock_tube import ShockTube
from PyShockflow.config import Config

config = Config('input.ini')
tube = ShockTube(config)
tube.SolveSystem()
