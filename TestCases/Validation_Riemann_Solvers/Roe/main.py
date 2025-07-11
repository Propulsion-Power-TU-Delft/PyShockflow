from PyShockflow.shock_tube import ShockTube
from PyShockflow.config import Config

testNumbers = [1,3,4,5]
inputFiles = ['input_Test%i.ini' % i for i in testNumbers]
for inputFile in inputFiles:
    config = Config(inputFile)
    tube = ShockTube(config)
    tube.SolveSystem()
