from pyshockflow import Driver
import numpy as np
import matplotlib.pyplot as plt
import pickle
# from PyShockflow.styles import *
from scipy.optimize import fsolve
from scipy.optimize import brentq

# INPUT
inputPkls = [
                'Results/test_coolprop_NX_500/Results.pik',
                'Results/test_refprop_NX_500/Results.pik',
                'Results/test_stanmix_NX_500/Results.pik',
            ]

soundSpeed = 34.21  # sound speed in the left state
rswSpeed = 35       # expected speed of the rarefaction shock wave

labels = [
            'CoolProp', 
            'RefProp',
            'StanMix',
            ]

linestyles = ['-', '--', '-.']

datas = []
for inputPkl in inputPkls:    
    with open(inputPkl, 'rb') as file:
        data = pickle.load(file)
    datas.append(data)

time = datas[0]['Time']
lasttime = datas[0]['Time'][-1]

# plt.figure(figsize=(12, 8))
# for i in range(len(time)):
#     plt.clf()
#     for data, label, ls in zip(datas, labels, linestyles):
#         iTime = np.argmin(np.abs(data['Time']-time[i]))
#         plt.plot(data['X Coords'], data['Primitive']['Pressure'][:, iTime], ls, label=label)
#     plt.axvline(x=1-soundSpeed*time[i], color='k', linestyle='--', lw=0.75, label='Sonic Limit')
#     plt.axvline(x=1-rswSpeed*time[i], color='r', linestyle=':', lw=0.75, label='Expected RSW')
#     plt.legend(loc='upper right', ncols=2)
#     plt.xlabel('X [m]')
#     plt.ylabel('Pressure [bar]')
#     plt.pause(0.01)


# plt.figure(figsize=(12, 8))
# for i in range(len(time)):
#     plt.clf()
#     for data, label, ls in zip(datas, labels, linestyles):
#         iTime = np.argmin(np.abs(data['Time']-time[i]))
#         plt.plot(data['X Coords'], data['Primitive']['Density'][:, iTime], ls, label=label)
#     plt.axvline(x=1-soundSpeed*time[i], color='k', linestyle='--', lw=0.75, label='Sonic Limit')
#     plt.axvline(x=1-rswSpeed*time[i], color='r', linestyle=':', lw=0.75, label='Expected RSW')
#     plt.legend(loc='upper right', ncols=2)
#     plt.xlabel('X [m]')
#     plt.ylabel('Density [kg/m3]')
#     plt.pause(0.01)
    

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
(ax1, ax2, ax3) = axes

# --- Density ---
for data, label, ls in zip(datas, labels, linestyles):
    ax1.plot(data['X Coords'], data['Primitive']['Density'][:, -1], ls, label=label)
ax1.axvline(x=1 - soundSpeed*lasttime, color='k', linestyle='--', lw=1, label='Sonic Limit')
ax1.axvline(x=1 - rswSpeed*lasttime, color='r', linestyle=':', lw=1, label='Expected RSW')
ax1.set_xlabel(r'$x$ [m]')
ax1.set_ylabel(r'$\rho$ [kg/m$^3$]')
fig.legend(loc='upper center', ncol=len(labels)+2, bbox_to_anchor=(0.53, 1.1))

# --- Velocity (middle) ---
for data, label, ls in zip(datas, labels, linestyles):
    ax2.plot(data['X Coords'], data['Primitive']['Velocity'][:, -1], ls, label=label)
ax2.axvline(x=1 - soundSpeed*lasttime, color='k', linestyle='--', lw=1)
ax2.axvline(x=1 - rswSpeed*lasttime, color='r', linestyle=':', lw=1)
ax2.set_xlabel(r'$x$ [m]')
ax2.set_ylabel(r'$u$ [m/s]')

# --- Pressure (last) ---
for data, label, ls in zip(datas, labels, linestyles):
    ax3.plot(data['X Coords'], data['Primitive']['Pressure'][:, -1] / 1e5, ls, label=label)
ax3.axvline(x=1 - soundSpeed*lasttime, color='k', linestyle='--', lw=1)
ax3.axvline(x=1 - rswSpeed*lasttime, color='r', linestyle=':', lw=1)
ax3.set_xlabel(r'$x$ [m]')
ax3.set_ylabel(r'$p$ [bar]')

for ax in axes:
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 1)

fig.tight_layout()
plt.savefig('FAST_RSW_Allresults.pdf', bbox_inches='tight')



plt.show()