import matplotlib.pyplot as plt
import numpy as np
import pickle
from NumericalCodes.riemann_problem import RiemannProblem
from NumericalCodes.shock_tube import ShockTube
import os


refFiles = ["../RiemannProblem_Analytical/solutions/Test%i.pik" %i for i in range(1,6)]
godFiles = ["solutions_god/Test%i.pik" %i for i in range(1,6)]
wafFiles = ["solutions_WAF/Test%i.pik" %i for i in range(1,6)]
roeFiles = ["solutions_Roe/Test%i.pik" %i for i in range(1,6)]


outFolder = 'pictures'
os.makedirs(outFolder, exist_ok=True)

for i in [0]:
    
    # reference file
    with open(refFiles[i], 'rb') as file:
        ref = pickle.load(file)
    
    # godunov file
    with open(godFiles[i], 'rb') as file:
        god = pickle.load(file)
    
    with open(wafFiles[i], 'rb') as file:
        waf = pickle.load(file)
    
    with open(roeFiles[i], 'rb') as file:
        roe = pickle.load(file)
    
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    ax[0, 0].plot(ref.x+0.5, ref.rho[:, -1], '-C0o', ms=2, label='Reference') # the analytical solution needs to be shifted by 0.5 in x
    ax[0, 0].plot(god.xNodes, god.solution['Density'][1:-1, -1], '-C1o', ms=2, mfc='none', label='Godunov', lw=0.75)  # don't plot the halo nodes
    # ax[0, 0].plot(waf.xNodes, waf.solution['Density'][1:-1, -1], '-C2o', ms=2, mfc='none', label='WAF', lw=0.75)  # don't plot the halo nodes
    ax[0, 0].plot(roe.xNodes, roe.solution['Density'][1:-1, -1], '-C3o', ms=2, mfc='none', label='Roe', lw=0.75)  # don't plot the halo nodes
    ax[0, 0].set_ylabel(r'Density')

    ax[0, 1].plot(ref.x+0.5, ref.u[:, -1], '-C0o', ms=2, label='Reference')
    ax[0, 1].plot(god.xNodes, god.solution['Velocity'][1:-1, -1], '-C1o', ms=2, mfc='none', label='Godunov', lw=0.75)
    # ax[0, 1].plot(waf.xNodes, waf.solution['Velocity'][1:-1, -1], '-C2o', ms=2, mfc='none', label='WAF', lw=0.75)
    ax[0, 1].plot(roe.xNodes, roe.solution['Velocity'][1:-1, -1], '-C3o', ms=2, mfc='none', label='Roe', lw=0.75)
    ax[0, 1].set_ylabel(r'Velocity')

    ax[1, 0].plot(ref.x+0.5, ref.p[:, -1], '-C0o', ms=2, label='Reference')
    ax[1, 0].plot(god.xNodes, god.solution['Pressure'][1:-1, -1], '-C1o', ms=2, mfc='none', label='Godunov', lw=0.75)
    # ax[1, 0].plot(waf.xNodes, waf.solution['Pressure'][1:-1, -1], '-C2o', ms=2, mfc='none', label='WAF', lw=0.75)
    ax[1, 0].plot(roe.xNodes, roe.solution['Pressure'][1:-1, -1], '-C3o', ms=2, mfc='none', label='Roe', lw=0.75)
    ax[1, 0].set_ylabel(r'Pressure')

    ax[1, 1].plot(ref.x+0.5, ref.e[:, -1], '-C0o', ms=2, label='Reference')
    ax[1, 1].plot(god.xNodes, god.solution['Energy'][1:-1, -1], '-C1o', ms=2, mfc='none', label='Godunov', lw=0.75)
    # ax[1, 1].plot(waf.xNodes, waf.solution['Energy'][1:-1, -1], '-C2o', ms=2, mfc='none', label='WAF', lw=0.75)
    ax[1, 1].plot(roe.xNodes, roe.solution['Energy'][1:-1, -1], '-C3o', ms=2, mfc='none', label='Roe', lw=0.75)
    ax[1, 1].set_ylabel(r'Energy')

    fig.suptitle('Test %i, Time %.3f' % (i+1, roe.timeVec[-1]))

    for row in ax:
            for col in row:
                col.set_xlabel('x')
                col.grid(alpha=.3)
                col.legend()
    
    plt.savefig(outFolder + '/Test%i.pdf' % (i+1), bbox_inches='tight')



plt.show()