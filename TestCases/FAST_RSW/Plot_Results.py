from PyShockflow.shock_tube import ShockTube
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
            'CoolProp-500', 
            'RefProp-500',
            'StanMix-500',
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
    

plt.figure()
for data, label, ls in zip(datas, labels, linestyles):
    plt.plot(data['X Coords'], data['Primitive']['Density'][:, -1], ls, label=label)
plt.axvline(x=1-soundSpeed*lasttime, color='k', linestyle='--', lw=0.75, label='Sonic Limit')
plt.axvline(x=1-rswSpeed*lasttime, color='r', linestyle=':', lw=0.75, label='Expected RSW')
plt.legend()
plt.xlabel('X [m]')
plt.ylabel('Density [kg/m3]')
plt.savefig('Density.pdf', bbox_inches='tight')

plt.figure()
for data, label, ls in zip(datas, labels, linestyles):
    plt.plot(data['X Coords'], data['Primitive']['Pressure'][:, -1]/1e5, ls, label=label)
plt.axvline(x=1-soundSpeed*lasttime, color='k', linestyle='--', lw=0.75, label='Sonic Limit')
plt.axvline(x=1-rswSpeed*lasttime, color='r', linestyle=':', lw=0.75, label='Expected RSW')
plt.legend()
plt.xlabel('X [m]')
plt.ylabel('Pressure [bar]')
plt.savefig('Pressure.pdf', bbox_inches='tight')

plt.figure()
for data, label, ls in zip(datas, labels, linestyles):
    plt.plot(data['X Coords'], data['Primitive']['Velocity'][:, -1], ls, label=label)
plt.axvline(x=1-soundSpeed*lasttime, color='k', linestyle='--', lw=0.75, label='Sonic Limit')
plt.axvline(x=1-rswSpeed*lasttime, color='r', linestyle=':', lw=0.75, label='Expected RSW')
plt.legend()
plt.xlabel('X [m]')
plt.ylabel('Velocity [m/s]')
plt.savefig('Velocity.pdf', bbox_inches='tight')


plt.show()