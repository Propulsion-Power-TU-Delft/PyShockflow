from pyshockflow import ShockTube
from pyshockflow import Config

config = Config('input_real_vinokur.ini')
tube = ShockTube(config)
tube.SolveSystem()
