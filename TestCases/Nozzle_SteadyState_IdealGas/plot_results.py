import numpy as np
import matplotlib.pyplot as plt
import pickle
from PyShockTube.styles import *

pressureList = [45, 75, 90, 94, 97]
pickleList = ['Results/outletPressure_%ikPa_NX_200/Results.pik' %pressure for pressure in pressureList]


fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()
fig4, ax4 = plt.subplots()
fig5, ax5 = plt.subplots()

for i, pickleFile in enumerate(pickleList):
    with open(pickleFile, 'rb') as file:
        solution = pickle.load(file)
    
    xCoords = solution['X Coords'][1:-1]
    density = solution['Primitive']["Density"][1:-1,-1]
    pressure = solution['Primitive']["Pressure"][1:-1,-1]
    velocity = solution['Primitive']["Velocity"][1:-1,-1]
    mach = solution['Fluid'].ComputeMach_u_p_rho(velocity, pressure, density)
    entropy = solution['Fluid'].ComputeEntropy_p_rho(pressure, density)
    totalPressure = solution['Fluid'].ComputeTotalPressure_p_M(pressure, mach)
    temperature = solution['Fluid'].ComputeTemperature_p_rho(pressure, density)
    totalTemperature = solution['Fluid'].ComputeTotalTemperature_T_M(temperature, mach)
    
    
    ax1.plot(xCoords, mach, label=r'$p_{out}=%i$ kPa' %(pressureList[i]))
    ax1.set_ylabel(r'Mach [-]')
    
    ax2.plot(xCoords, pressure/1e3, label=r'$p_{out}=%i$ kPa' %(pressureList[i]))
    ax2.set_ylabel(r'Pressure [kPa]')
    
    ax3.plot(xCoords, entropy/1e3, label=r'$p_{out}=%i$ kPa' %(pressureList[i]))
    ax3.set_ylabel(r'Entropy [kJ/kgK]')
    
    ax4.plot(xCoords, totalPressure/1e3, label=r'$p_{out}=%i$ kPa' %(pressureList[i]))
    ax4.set_ylabel(r'Total Pressure [kPa]')
    
    ax5.plot(xCoords, totalTemperature, label=r'$p_{out}=%i$ kPa' %(pressureList[i]))
    ax5.set_ylabel(r'Total Temperature [K]')
    
    for ax in [ax1, ax2, ax3, ax4, ax5]:
        ax.set_xlabel(r'$x$ [-]')
        ax.legend()
        ax.grid(alpha=.3)

fig1.savefig('Pictures/mach.pdf', bbox_inches='tight')
fig2.savefig('Pictures/pressure.pdf', bbox_inches='tight')
    
plt.show()
        
        
        
    