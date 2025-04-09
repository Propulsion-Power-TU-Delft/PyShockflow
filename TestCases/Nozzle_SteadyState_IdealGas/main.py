from PyShockTube.shock_tube import ShockTube
from PyShockTube.config import Config

pressureList = [97]
configList = ['input_%ikPa.ini' %(pressure) for pressure in pressureList]

for configFile in configList:
    config = Config(configFile)
    tube = ShockTube(config)
    tube.SolveSystem()
