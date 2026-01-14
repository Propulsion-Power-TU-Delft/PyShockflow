from pyshockflow import ShockTube
from pyshockflow import Config

inputFiles = ['input_Test%i.ini' % i for i in range(3, 5)]
for inputFile in inputFiles:
    config = Config(inputFile)
    tube = ShockTube(config)
    tube.SolveSystem()
