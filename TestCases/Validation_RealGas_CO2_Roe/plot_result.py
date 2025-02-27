import matplotlib.pyplot as plt
import numpy as np
import pickle
from PyShockTube.riemann_problem import RiemannProblem
from PyShockTube.shock_tube import ShockTube
import os

solFile = "solutions/co2_noentropyfix_500.pik"
solFile_fix = "solutions/co2_entropyfix_NX_100_TMAX_0.050000.pik"
solFile_high = "solutions/co2_entropyfix_100_highOrder.pik"
outFolder = 'Pictures'
os.makedirs(outFolder, exist_ok=True)

rhoRef = np.loadtxt("ReferenceData/rho.txt", skiprows=1)
pRef = np.loadtxt("ReferenceData/p.txt", skiprows=1)
uRef = np.loadtxt("ReferenceData/u.txt", skiprows=1)

# reference file
with open(solFile, 'rb') as file:
    sol = pickle.load(file)
with open(solFile_fix, 'rb') as file:
    solFix = pickle.load(file)
with open(solFile_high, 'rb') as file:
    solHigh = pickle.load(file)

fig, ax = plt.subplots(1, 3, figsize=(17, 6))

ax[0].plot(sol.xNodes, sol.solution['Density'][1:-1, -1], '-C0', label='Current') 
ax[0].plot(solFix.xNodes, solFix.solution['Density'][1:-1, -1], '-C1', ms=1.5, label=r'Current - $s$-fix') 
ax[0].plot(rhoRef[:,0], rhoRef[:,1], 'C2o', mfc='none', ms=5, label="D'Alessandro et al.")
# ax[0].plot(solHigh.xNodes, solHigh.solution['Density'][1:-1, -1], '-C4', ms=1.5, label=r'Current - $s$-fix - MUSCL') 
ax[0].set_title(r'$\rho \ \rm{[kg/m^3]}$')

ax[1].plot(sol.xNodes, sol.solution['Pressure'][1:-1, -1]/1e6, '-C0')
ax[1].plot(solFix.xNodes, solFix.solution['Pressure'][1:-1, -1]/1e6, '-C1')
ax[1].plot(pRef[:,0], pRef[:,1], 'C2o', mfc='none', ms=5)
# ax[1].plot(solHigh.xNodes, solHigh.solution['Pressure'][1:-1, -1]/1e6, '-C4', ms=1.5)
ax[1].set_title(r'$P \ \rm{[MPa]}$')

ax[2].plot(sol.xNodes, sol.solution['Velocity'][1:-1, -1], '-C0') 
ax[2].plot(solFix.xNodes, solFix.solution['Velocity'][1:-1, -1], '-C1')
ax[2].plot(uRef[:,0], uRef[:,1], 'C2o', mfc='none', ms=5) 
# ax[2].plot(solHigh.xNodes, solHigh.solution['Velocity'][1:-1, -1], '-C4', ms=1.50)
ax[2].set_title(r'$u \ \rm{[m/s]}$')

for row in ax:
    row.set_xlabel(r'$x \ \rm{[m]}$')
    row.grid(alpha=.3)
fig.legend(loc='lower center', bbox_to_anchor=(0.5, -0.05), ncol=3)

plt.savefig(outFolder + '/co2.pdf', bbox_inches='tight')



plt.show()