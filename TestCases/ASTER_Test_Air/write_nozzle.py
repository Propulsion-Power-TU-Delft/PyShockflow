import numpy as np
import matplotlib.pyplot as plt

"""
GENERATE CSV FILE FOR A NOZZLE DESCRIBED BY A PARABOLIC SHAPE, WITH INLET AND OUTLET SECTION AREA OF SIZE 1, AND THROTTLE SECTION IN THE MIDDLE BETWEEN INLET AND OUTLET
"""

xStart = 2.95            # axial coordinate of nozzle inlet (same reference of the tube)
xEnd = 3.05              # axial coordinate of nozzle outlet (same reference of the tube)
A_throat = 0.5          # throat area of the nozzle (normalized by area of the tube)
nPoints = 100            # number of points to describe the nozzle shape in the csv file


x = np.linspace(xStart, xEnd, nPoints)
Acoeff = -4*(A_throat-1)
Bcoeff = -Acoeff
Ccoeff = 1
z = (x-x[0])/(x[-1]-x[0])
Area = Acoeff*z**2 + Bcoeff*z + Ccoeff

# PLOT THE NOZZLE SHAPE
plt.figure()
plt.plot(x, Area)
plt.xlabel('x [m]')
plt.ylabel(r'$A/A_{tube}$ [-]')
plt.grid(alpha=.3)

with open('nozzle_%i.csv' %(A_throat*100), 'w') as file:
    file.write('x,A\n')
    for i in range(len(x)):
        file.write('%.6f,%.6f\n' %(x[i], Area[i]))

plt.show()



