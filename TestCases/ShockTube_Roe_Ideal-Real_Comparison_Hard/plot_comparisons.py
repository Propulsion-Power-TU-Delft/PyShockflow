import matplotlib.pyplot as plt
import numpy as np
import pickle
from NumericalCodes.riemann_problem import RiemannProblem
from NumericalCodes.shock_tube import ShockTube
import os


idealFile = "solutions/ideal_gas_0.015.pik"
realFile = "solutions/real_gas_0.015.pik"
outFolder = 'Pictures'
os.makedirs(outFolder, exist_ok=True)


# reference file
with open(idealFile, 'rb') as file:
    ideal = pickle.load(file)

# godunov file
with open(realFile, 'rb') as file:
    real = pickle.load(file)

fig, ax = plt.subplots(2, 2, figsize=(12, 8))
ax[0, 0].plot(ideal.xNodes, ideal.solution['Density'][1:-1, -1], '-C0o', ms=2, label='ideal gas') # the analytical solution needs to be shifted by 0.5 in x
ax[0, 0].plot(real.xNodes, real.solution['Density'][1:-1, -1], '-C1o', ms=2, mfc='none', label='real gas')  # don't plot the halo nodes
ax[0, 0].set_ylabel(r'Density')

ax[0, 1].plot(ideal.xNodes, ideal.solution['Velocity'][1:-1, -1], '-C0o', ms=2, label='ideal gas')
ax[0, 1].plot(real.xNodes, real.solution['Velocity'][1:-1, -1], '-C1o', ms=2, mfc='none', label='real gas')
ax[0, 1].set_ylabel(r'Velocity')

ax[1, 0].plot(ideal.xNodes, ideal.solution['Pressure'][1:-1, -1], '-C0o', ms=2, label='ideal gas')
ax[1, 0].plot(real.xNodes, real.solution['Pressure'][1:-1, -1], '-C1o', ms=2, mfc='none', label='real gas')
ax[1, 0].set_ylabel(r'Pressure')

ax[1, 1].plot(ideal.xNodes, ideal.solution['Energy'][1:-1, -1], '-C0o', ms=2, label='ideal gas')
ax[1, 1].plot(real.xNodes, real.solution['Energy'][1:-1, -1], '-C1o', ms=2, mfc='none', label='real gas')
ax[1, 1].set_ylabel(r'Energy')

fig.suptitle('Time %.3f' % (ideal.timeVec[-1]))

for row in ax:
        for col in row:
            col.set_xlabel('x')
            col.grid(alpha=.3)
            col.legend()

plt.savefig(outFolder + '/Comparison.pdf', bbox_inches='tight')



plt.show()