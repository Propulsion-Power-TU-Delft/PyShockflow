import matplotlib.pyplot as plt
import numpy as np
import pickle
from PyShockTube.riemann_problem import RiemannProblem
from PyShockTube.shock_tube import ShockTube
import os

solFile = "solutions/WV1_NX_250_TMAX_0.000764.pik"
outFolder = 'Pictures'
os.makedirs(outFolder, exist_ok=True)

rhoRef = np.loadtxt("ReferenceData/Test_WV1/rho.txt", skiprows=1)
pRef = np.loadtxt("ReferenceData/Test_WV1/p.txt", skiprows=1)

LENGTH = 1
NX = 150
UL, UR = 0.0, 0.0
Rhocr = 322
RHOL = 1.01*Rhocr
RHOR = 0.59400*Rhocr
Pcr = 22.064e6
PL, PR = 1.6077*Pcr, 0.8957*Pcr
TIME_MAX = 0.2*LENGTH*(Rhocr/Pcr)**(0.5)



# reference file
with open(solFile, 'rb') as file:
    sol = pickle.load(file)

fig, ax = plt.subplots(1, 2, figsize=(12, 6))

ax[0].plot(sol.xNodes, sol.solution['Density'][1:-1, -1]/Rhocr, '-C0', ms=1.5, label='Current') 
ax[0].plot(rhoRef[:,0],rhoRef[:,1], 'C1o', mfc='none', ms=5, label='Guardone et al.') 
ax[0].set_title(r'$\rho / \rho_{cr}$ [-]')

ax[1].plot(sol.xNodes, sol.solution['Pressure'][1:-1, -1]/Pcr, '-C0', label='Current') 
ax[1].plot(pRef[:,0],pRef[:,1], 'C1o', mfc='none', ms=5, label='Guardone et al.') 
ax[1].set_title(r'$P/P_{cr}$ [-]')

for row in ax:
    row.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    row.set_xlim([0,1])
    row.set_xlabel('x')
    row.grid(alpha=.3)
    row.legend()

plt.suptitle("Test WV1")
plt.savefig(outFolder + '/WV1.pdf', bbox_inches='tight')



plt.show()