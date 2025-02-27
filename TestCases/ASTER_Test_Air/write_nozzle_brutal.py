import numpy as np
import matplotlib.pyplot as plt

"""
GENERATE CSV FILE FOR A NOZZLE RESEMBLING THE THROAT OF ASTER SHAPE. THE AREA VALUES ARE NORMALIZED BY THE TUBE AREA, SO THAT IN THE SIMULATION REFERENCE_AREA CAN BE LEFT TO ONE.
"""

xThroat = 3.005                           # throat axial location
xStart = xThroat - (7.188)*1E-3     # axial coordinate of nozzle inlet 
xEnd = xThroat + 23.412*1E-3            # axial coordinate of nozzle outlet (same reference of the tube)
A_throat = (22/30.3)**2                 # throat area of the nozzle (normalized by area of the tube)
A_start = A_throat                           # throat area of the nozzle inlet (normalized by area of the tube)
A_end = A_throat                             # throat area of the nozzle outlet (normalized by area of the tube)
nPoints = 250                           # number of points to describe the nozzle shape in the csv file





x = np.linspace(xStart, xEnd, nPoints)
Area = np.interp(x, [xStart, xThroat, xEnd], [A_start, A_throat, A_end])

# PLOT THE NOZZLE SHAPE
plt.figure()
plt.plot([xStart, xThroat, xEnd], [A_start, A_throat, A_end], 'o', label='Control Points')
plt.plot(x, Area, label='Linear Interpolation')
plt.xlabel('x [m]')
plt.ylabel(r'$A/A_{tube}$ [-]')
plt.grid(alpha=.3)
plt.legend()

with open('nozzle.csv', 'w') as file:
    file.write('x,A\n')
    for i in range(len(x)):
        file.write('%.6f,%.6f\n' %(x[i], Area[i]))

plt.show()



