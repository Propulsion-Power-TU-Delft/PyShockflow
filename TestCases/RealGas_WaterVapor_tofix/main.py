from PyShockTube.shock_tube import ShockTube
from PyShockTube.config import Config

config = Config('input_DG1.ini')
tube = ShockTube(config)
tube.SolveSystem()
