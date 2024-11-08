import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
sys.path.append("../")
from src.riemann_problem import RiemannProblem
from src.shock_tube import ShockTube


filepath = 'solutions/Test2.pik'
with open(filepath, 'rb') as file:
    tube = pickle.load(file)

tube.ShowAnimation()
tube.PlotSolution(iTime=tube.nTime-1)



plt.show()