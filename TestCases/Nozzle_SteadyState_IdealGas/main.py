from PyShockflow.shock_tube import ShockTube
from PyShockflow.config import Config

pressureList = [45, 75, 90, 94]
configList = ['input_%ikPa.ini' %(pressure) for pressure in pressureList]

for configFile in configList:
    config = Config(configFile)
    tube = ShockTube(config)
    tube.SolveSystem()
