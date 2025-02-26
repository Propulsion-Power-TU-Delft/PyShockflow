from PyShockTube.shock_tube import ShockTube
from PyShockTube.config import Config
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

with open('Results/ASTER_AIR_NX_250_TMAX_0.005000.pik', 'rb') as file:
    tube = pickle.load(file)
    
gmma = 1.4
At = 0.7
Rgas = 287.14

axialCoord = tube.xNodes
time = tube.timeVec
pressure = tube.solution['Pressure'][1:-1,:]
density = tube.solution['Density'][1:-1,:]
temperature = pressure/density/Rgas
velocity = tube.solution['Velocity'][1:-1,:]
mach = velocity/np.sqrt(gmma*pressure/density)

ni,nt = mach.shape
id1 = int(0.8/2*ni) # 1.2m from throat located at 2 m
id2 = int(0.6/2*ni) # 1.6m from throat located at 2 m
id3 = int(0.4/2*ni) # 1.8m from throat located at 2 m
id4 = int(0.2/2*ni) # 1.8m from throat located at 2 m

plt.figure()
plt.plot(time, pressure[id1,:]/1e5, label='kulite 1')
plt.plot(time, pressure[id2,:]/1e5, label='kulite 2')
plt.plot(time, pressure[id3,:]/1e5, label='kulite 3')
plt.plot(time, pressure[id4,:]/1e5, label='kulite 3')
plt.xlabel(r'$t \ \rm{[s]}$')
plt.ylabel(r'$p \ \rm{[bar]}$')
plt.grid(alpha=.3)
plt.legend()






# CHECK THE MACH NUMBER
machThroat = mach[np.argmin(tube.areaTube), :]
areaRatioExit = tube.areaTube[-1]/np.min(tube.areaTube)
def compute_alpha_residual(M):
    return areaRatioExit-1/M*(2/(gmma+1)*(1+(gmma-1)/2*M**2))**((gmma+1)/(2*gmma-2))

machExitTheoretical = fsolve(compute_alpha_residual, 1.5)
plt.figure()
plt.plot(time*1e3, machThroat, 'C0',label='Throat Section')
plt.plot(time*1e3, np.zeros_like(time)+1, '--C0', label=r'Theoretical Steady Value: $M=1$')
plt.plot(time*1e3, mach[-1,:], 'C1', label='Exit Section')
plt.plot(time*1e3, np.zeros_like(time) + machExitTheoretical, '--C1', label=r'Theoretical Steady Value: $\alpha(M) = \frac{A}{A^*}$')
plt.xlabel(r'$t \ \rm{[ms]}$')
plt.ylabel(r'$M \ \rm{[-]}$')
plt.legend()
plt.grid(alpha=.3)
plt.savefig('mach.pdf', bbox_inches='tight')


# CHECK THE MASS FLOW RATE
massflux = density[-1,:]*velocity[-1,:] # numerical value
totPressure = pressure*(1+(gmma-1)/2*mach**2)**((gmma)/(gmma-1))
totTemperature = temperature*(1+(gmma-1)/2*mach**2)
massfluxTheoretical = totPressure*np.min(tube.areaTube)/(np.sqrt(gmma*Rgas*totTemperature))*gmma*(2/(gmma+1))**((gmma+1)/(2*gmma-2))

plt.figure()
plt.plot(time*1e3, massflux, label='Simulation')
plt.plot(time*1e3, np.zeros_like(time) + massfluxTheoretical[-1,-1], '--r', label=r'Theoretical Steady Value: $\dot{m}=\frac{P_t A^*}{\sqrt{\gamma R T_t}} \cdot f(\gamma)$')
plt.xlabel(r'$t \ \rm{[ms]}$')
plt.ylabel(r'$\dot{m}_{EXIT} \ \rm{[kg/s]}$')
plt.legend()
plt.grid(alpha=.3)
plt.savefig('massflow.pdf', bbox_inches='tight')

plt.show()
