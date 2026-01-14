from pyshockflow import ShockTube
from pyshockflow import Config

config = Config('input.ini')
tube = ShockTube(config)
tube.SolveSystem()
