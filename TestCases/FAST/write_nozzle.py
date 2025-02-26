import numpy as np
import matplotlib.pyplot as plt

xStart = 1.0 # assuming a tube long 5 meters for the moment
xEnd = 1.0+0.15/2+0.5 # half nozzle plus 0.5 meters of variation to tank area
A_tube = 1256/1256
A_throat = 531.2/1256 # ratio between the throat area of the nozzle and tube normal area
A_tank = 1E6/1256 # 1 m2 for tank area size
x_tube = xStart
x_throat = xStart + 0.075
x_tank = xEnd

nPoints = 100 # to describe the nozzle


# use polyfit polyval
x = np.linspace(xStart, xEnd, nPoints)
Area = np.interp(x, [x_tube, x_throat, x_tank], [A_tube, A_throat, A_tank])

plt.figure()
plt.plot(x, Area)
plt.xlabel('x')
plt.ylabel('A')

with open('nozzle.csv', 'w') as file:
    file.write('x,A\n')
    for i in range(len(x)):
        file.write('%.6f,%.6f\n' %(x[i], Area[i]))

plt.show()



