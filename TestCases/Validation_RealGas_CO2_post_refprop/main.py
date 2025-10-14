from PyShockflow.shock_tube import ShockTube
from PyShockflow.config import Config
# import sys

# sys.path.append('/Users/fneri/Documents/PhD/fluidproperties')
# import fluid_properties

config = Config('input_real_arabi.ini')
tube = ShockTube(config)
tube.SolveSystem()
