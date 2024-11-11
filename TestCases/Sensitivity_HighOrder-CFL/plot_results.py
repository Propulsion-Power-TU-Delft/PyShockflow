import matplotlib.pyplot as plt
import numpy as np
import pickle
from PyShockTube.riemann_problem import RiemannProblem
from PyShockTube.shock_tube import ShockTube
import os


refFile = "../IdealGas_Analytical/solutions/Test1.pik"
cfls = [0.01, 0.05, 0.1, 0.20]
roeFiles = ["solutions/SodsTest_%.2f.pik" %i for i in cfls]
outFolder = 'pictures'
os.makedirs(outFolder, exist_ok=True)


# reference file
with open(refFile, 'rb') as file:
    ref = pickle.load(file)

roes = []
for input_file in roeFiles:
    with open(input_file, 'rb') as file:
        roes.append(pickle.load(file))


fig, ax = plt.subplots(2, 2, figsize=(12, 8))
ax[0, 0].plot(ref.x+0.5, ref.rho[:, -1], '-C0', ms=2, label='Reference') # the analytical solution needs to be shifted by 0.5 in x
for i, roe in enumerate(roes):
    ax[0, 0].plot(roe.xNodes, roe.solution['Density'][1:-1, -1], '--C%i' %(i+1), mfc='none', label='CFL: %.2f' %cfls[i])  # don't plot the halo nodes
ax[0, 0].set_ylabel(r'Density')

ax[0, 1].plot(ref.x+0.5, ref.u[:, -1], '-C0', ms=2, label='Reference')
for i, roe in enumerate(roes):
    ax[0, 1].plot(roe.xNodes, roe.solution['Velocity'][1:-1, -1], '--C%i' %(i+1), mfc='none', label='CFL: %.2f' %cfls[i])
ax[0, 1].set_ylabel(r'Velocity')

ax[1, 0].plot(ref.x+0.5, ref.p[:, -1], '-C0', ms=2, label='Reference')
for i, roe in enumerate(roes):
    ax[1, 0].plot(roe.xNodes, roe.solution['Pressure'][1:-1, -1], '--C%i' %(i+1), mfc='none', label='CFL: %.2f' %cfls[i])
ax[1, 0].set_ylabel(r'Pressure')

ax[1, 1].plot(ref.x+0.5, ref.e[:, -1], '-C0', ms=2, label='Reference')
for i, roe in enumerate(roes):
    ax[1, 1].plot(roe.xNodes, roe.solution['Energy'][1:-1, -1], '--C%i' %(i+1), mfc='none', label='CFL: %.2f' %cfls[i])
ax[1, 1].set_ylabel(r'Energy')


for row in ax:
        for col in row:
            col.set_xlabel('x')
            col.grid(alpha=.3)
            col.legend()

plt.savefig(outFolder + '/CFL_Sensitivity.pdf', bbox_inches='tight')



plt.show()