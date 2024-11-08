import matplotlib.pyplot as plt
import numpy as np
import pickle
from NumericalCodes.riemann_problem import RiemannProblem
from NumericalCodes.shock_tube import ShockTube
import os


idealFile = "solutions/ideal_0.002500.pik"
realFile = "solutions/real_0.002500.pik"
outFolder = 'Pictures'
os.makedirs(outFolder, exist_ok=True)


# reference file
with open(idealFile, 'rb') as file:
    ideal = pickle.load(file)

# godunov file
with open(realFile, 'rb') as file:
    real = pickle.load(file)

fig, ax = plt.subplots(2, 2, figsize=(12, 8))

for iTime in range(ideal.timeVec.shape[0]):

    for row in ax:
            for col in row:
                col.cla()

    ax[0, 0].plot(ideal.xNodes, ideal.solution['Density'][1:-1, iTime], '-C0o', ms=2, label='ideal gas') 
    ax[0, 0].plot(real.xNodes, real.solution['Density'][1:-1, iTime], '-C1o', ms=2, mfc='none', label='real gas')  
    ax[0, 0].set_ylabel(r'Density [kg/m3]')

    ax[0, 1].plot(ideal.xNodes, ideal.solution['Velocity'][1:-1, iTime], '-C0o', ms=2, label='ideal gas')
    ax[0, 1].plot(real.xNodes, real.solution['Velocity'][1:-1, iTime], '-C1o', ms=2, mfc='none', label='real gas')
    ax[0, 1].set_ylabel(r'Velocity [m/s]')

    ax[1, 0].plot(ideal.xNodes, ideal.solution['Pressure'][1:-1, iTime]/1e5, '-C0o', ms=2, label='ideal gas')
    ax[1, 0].plot(real.xNodes, real.solution['Pressure'][1:-1, iTime]/1e5, '-C1o', ms=2, mfc='none', label='real gas')
    ax[1, 0].set_ylabel(r'Pressure [bar]')

    ax[1, 1].plot(ideal.xNodes, (ideal.solution['Energy'][1:-1, iTime]-ideal.solution['Energy'][1, iTime])/1e3, '-C0o', ms=2, label='ideal gas')
    ax[1, 1].plot(real.xNodes, (real.solution['Energy'][1:-1, iTime]-real.solution['Energy'][1, iTime])/1e3, '-C1o', ms=2, mfc='none', label='real gas')
    ax[1, 1].set_ylabel(r'Delta Energy [kJ/kg]')

    fig.suptitle('Time %.3f' % (ideal.timeVec[iTime]))

    for row in ax:
            for col in row:
                col.set_xlabel('x')
                col.grid(alpha=.3)
                col.legend()


    plt.pause(0.001)

# plt.savefig(outFolder + '/Comparison.pdf', bbox_inches='tight')



plt.show()