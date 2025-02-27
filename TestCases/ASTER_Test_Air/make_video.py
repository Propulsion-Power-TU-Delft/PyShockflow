import pickle 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os    
from PyShockTube.styles import *

#########################     INPUT         #########################

outputFolder = 'Videos'
picklePath = 'Results/ASTER_AIR_70_NOZZLE_REAL_NX_200_TMAX_0.050000.pik'            # Path to the pickle file of the simulation
maxFrames = 1000                            # choose how many snapshots you want to visualize (must be < than total snapshots of simulation)

#video settings
FPS = 30                                    # frames per second
DPI = 400                                   # definition (<500 works, more no, don't know why)

#####################################################################

os.makedirs(outputFolder, exist_ok=True)

# open the pickle file
with open(picklePath, 'rb') as file:
    solution = pickle.load(file)

# save aliases for the arrays of interest
x = solution.xNodes
time = solution.timeVec
density = solution.solution['Density'][1:-1,:]
velocity = solution.solution['Velocity'][1:-1,:]
pressure = solution.solution['Pressure'][1:-1,:]/1e5
staticEnergy = solution.solution['Energy'][1:-1,:]
mach = velocity/np.sqrt(1.4*pressure/density) # substitute this with real model calculation for real gas (tube.fluid.ComputeSoundSpeed.....)
nPoints, nTimes = density.shape
iterations = np.linspace(0, nTimes-1, num=maxFrames, dtype=int)

# fields = [rho, u, p, E, mach]
# labels = ['Density [kg/m3]', 'Velocity [m/s]', 'Pressure [Pa]', 'Static Energy [J/kg]' , 'Mach [-]']
# videoNames = ['Density.mp4', 'Velocity.mp4', 'Pressure.mp4', 'Energy.mp4', 'Mach.mp4']

fields = [mach, pressure]
labels = ['Mach [-]', 'Pressure [bar]']
videoNames = ['Mach.mp4', 'Pressure.mp4']

# PLOTS AND VIDEO
def plot_limits(f, extension=0.05):
    max = f.max()
    min = f.min()
    left = min-(max-min)*extension
    right = max+(max-min)*extension
    return left, right


for i,field in enumerate(fields):
    xmin, xmax = plot_limits(x)
    ymin, ymax = plot_limits(field)

    fig, ax = plt.subplots(constrained_layout=True)    
    line, = ax.plot([], [], '-C0')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r'$x$ [m]')
    ax.set_ylabel(labels[i])
    ax.grid(alpha=0.2)


    def update(iteration):
        line.set_data(x, field[:, iteration])
        ax.set_title(f'Time: {time[iteration]:.3e} [s]')
        return line,

    ani = animation.FuncAnimation(fig, update, frames=iterations, blit=False)

    # Save the animation as a video
    ani.save(outputFolder+"/"+videoNames[i], writer='ffmpeg', fps=FPS, dpi=DPI)
    print('Video %s Saved' %(videoNames[i]))

# plt.show()
