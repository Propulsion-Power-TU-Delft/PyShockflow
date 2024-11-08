import matplotlib.pyplot as plt
import numpy as np
import pickle
from NumericalCodes.riemann_problem import RiemannProblem
from NumericalCodes.shock_tube import ShockTube
import os


inputFile = "solutions/SodsTest_tMax_3.00.pik"

# reference file
with open(inputFile, 'rb') as file:
    res = pickle.load(file)

plt.figure()
for it in range(res.timeVec.shape[0]):
    plt.cla()
    plt.plot(res.xNodes, res.solution['Pressure'][1:-1, it])
    # plt.plot(res.xNodes+res.xNodes[-1], res.solution['Pressure'][1:-1, it])
    # plt.plot(res.xNodes+2*res.xNodes[-1], res.solution['Pressure'][1:-1, it])
    plt.xlabel('x')
    plt.ylabel('Pressure')
    plt.grid(alpha=0.3)
    plt.ylim([0, 1.05])
    plt.title('Time: %.3f' %res.timeVec[it])
    plt.pause(0.001)





plt.show()