from PyShockflow.shock_tube import ShockTube
from PyShockflow.config import Config

config = Config('input_real_arabi.ini')
tube = ShockTube(config)
tube.SolveSystem()
