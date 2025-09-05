from PyShockflow.shock_tube import ShockTube
from PyShockflow.config import Config

ins = ['input_normal.ini', 'input_refined.ini']

for input in ins:
    config = Config(input)
    tube = ShockTube(config)
    tube.SolveSystem()
