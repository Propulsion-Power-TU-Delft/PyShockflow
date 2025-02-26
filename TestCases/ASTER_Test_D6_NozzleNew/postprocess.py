from PyShockTube.shock_tube import ShockTube
from PyShockTube.config import Config
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

with open('Results/ASTER_D6_NX_250_TMAX_0.060000.pik', 'rb') as file:
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
id1 = int(0.67/1.8*ni) # 1.2m from diaphragm located at 1.87 m
id2 = int(0.27/1.8*ni) # 1.6m from diaphragm located at 1.87 m
id3 = int(0.07/1.8*ni) # 1.8m from diaphragm located at 1.87 m

plt.figure()
plt.plot(time, pressure[id1,:]/1e5, label='kulite 1')
plt.plot(time, pressure[id2,:]/1e5, label='kulite 2')
plt.plot(time, pressure[id3,:]/1e5, label='kulite 3')
plt.xlabel(r'$t \ \rm{[s]}$')
plt.ylabel(r'$p \ \rm{[bar]}$')
plt.grid(alpha=.3)
plt.legend()
plt.show()

