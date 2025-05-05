import matplotlib.pyplot as plt
import numpy as np
from PyShockflow.styles import *

residualsRho = []
residualsRhoU = []
residualsRhoE = []
with open('97kpa.log', 'r') as file:
    for line in file:
        if 'Residuals:' in line:
            tmp = line.split(':')[1].strip().split(',')
            tmp = [float(res) for res in tmp]
            residualsRho.append(tmp[0])
            residualsRhoU.append(tmp[1])
            residualsRhoE.append(tmp[2])

residualsRho = np.array(residualsRho)
residualsRhoU = np.array(residualsRhoU)
residualsRhoE = np.array(residualsRhoE)

def shift_to_drop(res):
    return res-res[0]

plt.figure()
plt.plot(shift_to_drop(residualsRho), label=r'$R(\rho)$')
plt.plot(shift_to_drop(residualsRhoU), label=r'$R(\rho u)$')
plt.plot(shift_to_drop(residualsRhoE), label=r'$R(\rho e_t)$')
plt.legend()
plt.grid(alpha=.3)
plt.xlabel('Iterations [-]')
plt.ylabel('Residuals Drop [-]')
plt.savefig('Pictures/residuals.pdf', bbox_inches='tight')



plt.show()