from PyShockflow.shock_tube import ShockTube
from PyShockflow.config import Config


nameList = ['1stOrder', 'minmod', 'vanalbada', 'vanleer', 'superbee']
inputList = ['input_%s.ini' %(name) for name in nameList]

for input in inputList:
    config = Config(input)
    tube = ShockTube(config)
    tube.SolveSystem()
