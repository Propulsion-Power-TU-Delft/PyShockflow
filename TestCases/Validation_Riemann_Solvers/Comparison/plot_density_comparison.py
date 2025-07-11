import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from PyShockflow.riemann_problem import RiemannProblem
from PyShockflow.styles import *

testNumber = [1,3,4,5]
markerSize = 4
# ANALYTICAL SOLUTIONS
analyticalResults = ['../Analytical_RiemannProblems/solutions/Test%i.pik' % i for i in testNumber]
godunovResults = ['../Godunov/Results/Test%i_NX_250/Results.pik' % i for i in testNumber]
roeResults = ['../Roe/Results/Test%i_NX_250/Results.pik' % i for i in testNumber]

figSize = (16, 4)
for iInput in range(len(analyticalResults)):
    
    fig, ax = plt.subplots(1,3, figsize=figSize)
    
    # ANALYTICAL
    with open(analyticalResults[iInput], 'rb') as file:
        res = pickle.load(file)
        ax[0].plot(res.x+0.5, res.rho[:,-1], label=r'Reference')
        ax[1].plot(res.x+0.5, res.u[:,-1])
        ax[2].plot(res.x+0.5, res.p[:,-1])
        
        ax[0].set_title(r'$\rho \ \rm{[-]}$')
        ax[1].set_title(r'$u \ \rm{[-]}$')
        ax[2].set_title(r'$p \ \rm{[-]}$')
        
    # GODUNOV
    with open(godunovResults[iInput], 'rb') as file:
        res = pickle.load(file)
        ax[0].plot(res['X Coords'][1:-1], res['Primitive']['Density'][1:-1,-1], 'o', ms=markerSize, mfc='none', label=r'Godunov')
        ax[1].plot(res['X Coords'][1:-1], res['Primitive']['Velocity'][1:-1,-1], 'o', ms=markerSize, mfc='none')
        ax[2].plot(res['X Coords'][1:-1], res['Primitive']['Pressure'][1:-1,-1], 'o', ms=markerSize, mfc='none')
    
    # ROE
    with open(roeResults[iInput], 'rb') as file:
        res = pickle.load(file)
        ax[0].plot(res['X Coords'][1:-1], res['Primitive']['Density'][1:-1,-1],'x', ms=markerSize, label=r'Roe')
        ax[1].plot(res['X Coords'][1:-1], res['Primitive']['Velocity'][1:-1,-1],'x', ms=markerSize)
        ax[2].plot(res['X Coords'][1:-1], res['Primitive']['Pressure'][1:-1,-1],'x', ms=markerSize)
        
    
    for axx in ax:
            axx.set_xlabel(r'$x \ \rm{[-]}$')
            axx.grid(alpha=.3)
    
    # Add a single legend at the bottom center
    # Add a single legend at the bottom center
    fig.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, +1.11))

    # Adjust layout to make room for the legend
    
    fig.savefig('Pictures/Test%i.pdf' % testNumber[iInput], bbox_inches='tight')


    
    

plt.show()
