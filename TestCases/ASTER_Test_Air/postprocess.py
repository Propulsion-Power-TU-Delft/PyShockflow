from PyShockTube.shock_tube import ShockTube
from PyShockTube.config import Config
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from PyShockTube.styles import *
import os

##########################################  INPUT ######################################################

pickleFile = 'Results/ASTER_AIR_75_NOZZLE_REAL_NX_200_TMAX_0.060000.pik'
outputFolder = 'Pictures'
suffixPictures = '75IdealGas' # suffix appended to pictures to distinguish different simulations

########################################################################################################


os.makedirs('Pictures', exist_ok=True)
with open(pickleFile, 'rb') as file:
    tube = pickle.load(file)


# ALIAS TO SOLUTION FIELDS (tube is the ShockTube object with the full space-time solution stored in it)
axialCoord = tube.xNodes
time = tube.timeVec
pressure = tube.solution['Pressure'][1:-1,:]    # avoid Halo nodes
density = tube.solution['Density'][1:-1,:]      # avoid Halo nodes
velocity = tube.solution['Velocity'][1:-1,:]    # avoid Halo nodes



# PLOT OF KULITE PRESSURE EVOLUTIONS IN TIME
ni,nt = velocity.shape
totalLength = 3.05              # total tube length (3m to nozzle throat, and other 5cm to exit in the tank)
id1 = int(1.8/totalLength*ni)   # 1.2m from throat located at 3m -> 1.8m from the left end
id2 = int(1.6/totalLength*ni)   # 1.4m from throat located at 3m -> 1.8m from the left end
id3 = int(1.4/totalLength*ni)   # 1.6m from throat located at 3m -> 1.8m from the left end
id4 = int(1.2/totalLength*ni)   # 1.8m from throat located at 3m -> 1.8m from the left end

plt.figure()
plt.plot(time*1e3, pressure[id1,:]/1e5, label='PT 1')
plt.plot(time*1e3, pressure[id2,:]/1e5, label='PT 2')
plt.plot(time*1e3, pressure[id3,:]/1e5, label='PT 3')
plt.plot(time*1e3, pressure[id4,:]/1e5, label='PT 3')
plt.xlabel(r'$t \ \rm{[ms]}$')
plt.ylabel(r'$p \ \rm{[bar]}$')
plt.grid(alpha=.3)
plt.legend()
plt.savefig(outputFolder+'/kulite_pressures_%s.pdf' %(suffixPictures), bbox_inches='tight')


# PLOT OF EXIT AND THROAT MACH NUMBERS IN TIME
machExit = np.zeros(nt)
machThroat = np.zeros(nt)
indexThroat = np.argmin(tube.areaTube)

for iTime in range(nt):
    machExit[iTime] = velocity[-1,iTime] / tube.fluid.ComputeSoundSpeed_p_rho(pressure[-1,iTime], density[-1, iTime])
    machThroat[iTime] = velocity[indexThroat,iTime] / tube.fluid.ComputeSoundSpeed_p_rho(pressure[indexThroat,iTime], density[indexThroat, iTime])

plt.figure()
plt.plot(time*1e3, machExit, label='Exit')
# plt.plot(time*1e3, machThroat, label='Throat')
plt.grid(alpha=grid_opacity)
plt.xlabel(r'$t \ \rm{[ms]}$')
plt.ylabel(r'$M$ [-]')
plt.legend()
plt.savefig(outputFolder+'/Mach_time_%s.pdf' %(suffixPictures), bbox_inches='tight')



# PLOT OF EXIT MASS FLUX IN TIME
massfluxExit = density[-1,:]*velocity[-1,:]
plt.figure()
plt.plot(time*1e3, massfluxExit, label='Exit')
plt.xlabel(r'$t \ \rm{[ms]}$')
plt.ylabel(r'$\frac{\dot{m}_{EXIT}}{A_{TUBE}} \ \rm{[kg/s/m^2]}$')
plt.legend()
plt.grid(alpha=.3)
plt.savefig(outputFolder+'/massflow_%s.pdf' %(suffixPictures), bbox_inches='tight')
plt.show()


