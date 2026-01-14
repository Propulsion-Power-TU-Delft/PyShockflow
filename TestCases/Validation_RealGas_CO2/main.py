from pyshockflow import ShockTube
from pyshockflow import Config

inputFiles = [
    'input_ideal_standard.ini',
    'input_ideal_vinokur.ini',
    'input_real_arabi.ini',
    'input_real_vinokur.ini'
]

for inputFile in inputFiles:
    config = Config(inputFile)
    tube = ShockTube(config)
    tube.solve()
