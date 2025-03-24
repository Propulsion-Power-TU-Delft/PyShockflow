import matplotlib.pyplot as plt
import numpy as np
import pickle
from PyShockTube.riemann_problem import RiemannProblem
from PyShockTube.shock_tube import ShockTube
import os


resultsFile = ["../IdealGas_Analytical/solutions/Test1.pik",
               "Results/NoRefinement_NX_250_TMAX_0.250000.pik",
               "Results/Refinement-04-06_NX_250_TMAX_0.250000.pik"]

labels = ['Reference', 'Normal Grid', 'Refined Grid']

outFolder = 'Pictures'
os.makedirs(outFolder, exist_ok=True)

resultsPickle = []
for i in range(len(resultsFile)):
    with open(resultsFile[i], 'rb') as file:
        resultsPickle.append(pickle.load(file))

reference = resultsPickle[0]
roeNormal = resultsPickle[1]
roeRefined = resultsPickle[2]

fig, ax = plt.subplots(2, 2, figsize=(12, 8))

ax[0, 0].plot(reference.x + 0.5, reference.rho[:, -1], '-C0', ms=2, label=labels[0])
ax[0, 0].plot(roeNormal.xNodesVirt, roeNormal.solution["Density"][:, -1], '-C1', ms=2, label=labels[1])
ax[0, 0].plot(roeRefined.xNodesVirt, roeRefined.solution["Density"][:, -1], '-C2', ms=2, label=labels[2])
ax[0, 0].set_ylabel(r'Density')

ax[0, 1].plot(reference.x + 0.5, reference.u[:, -1], '-C0', ms=2)
ax[0, 1].plot(roeNormal.xNodesVirt, roeNormal.solution["Velocity"][:, -1], '-C1')
ax[0, 1].plot(roeRefined.xNodesVirt, roeRefined.solution["Velocity"][:, -1], '-C2')
ax[0, 1].set_ylabel(r'Velocity')

ax[1, 0].plot(reference.x + 0.5, reference.p[:, -1], '-C0', ms=2)
ax[1, 0].plot(roeNormal.xNodesVirt, roeNormal.solution["Pressure"][:, -1], '-C1')
ax[1, 0].plot(roeRefined.xNodesVirt, roeRefined.solution["Pressure"][:, -1], '-C2')
ax[1, 0].set_ylabel(r'Pressure')

ax[1, 1].plot(reference.x + 0.5, reference.e[:, -1], '-C0', ms=2)
ax[1, 1].plot(roeNormal.xNodesVirt, roeNormal.solution["Energy"][:, -1], '-C1')
ax[1, 1].plot(roeRefined.xNodesVirt, roeRefined.solution["Energy"][:, -1], '-C2')
ax[1, 1].set_ylabel(r'Energy')

# Add legend only once for the figure
fig.legend(labels, loc='upper center', ncol=2)

for row in ax:
        for col in row:
            col.set_xlabel('x')
            col.grid(alpha=.2)

plt.savefig(outFolder + '/Comparison.pdf', bbox_inches='tight')



plt.show()