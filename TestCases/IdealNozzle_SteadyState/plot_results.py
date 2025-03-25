import numpy as np
import matplotlib.pyplot as plt
import pickle
from PyShockTube.styles import *

pressureList = [45,75,90,94,97]
pickleList = ['Results/outletPressure_%ikPa_NX_200_TMAX_0.050000.pik' %pressure for pressure in pressureList]


fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

for i, pickleFile in enumerate(pickleList):
    with open(pickleFile, 'rb') as file:
        solution = pickle.load(file)
    
    xCoords = solution.xNodes
    density = solution.solution["Density"][1:-1,-1]
    pressure = solution.solution["Pressure"][1:-1,-1]
    velocity = solution.solution["Velocity"][1:-1,-1]
    mach = solution.fluid.ComputeMach_u_p_rho(velocity, pressure, density)
    
    
    ax1.plot(xCoords, mach, label=r'$p_{out}=%i$ kPa' %(pressureList[i]))
    ax1.set_ylabel(r'Mach [-]')
    
    ax2.plot(xCoords, pressure/1e3, label=r'$p_{out}=%i$ kPa' %(pressureList[i]))
    ax2.set_ylabel(r'Pressure [kPa]')
    
    for ax in [ax1, ax2]:
        ax.set_xlabel(r'$x$ [-]')
        ax.legend()
        ax.grid(alpha=.3)

fig1.savefig('mach.pdf', bbox_inches='tight')
fig2.savefig('pressure.pdf', bbox_inches='tight')
    
plt.show()
        
        
        
    