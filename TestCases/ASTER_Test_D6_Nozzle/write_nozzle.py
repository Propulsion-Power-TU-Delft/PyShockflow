import numpy as np
import matplotlib.pyplot as plt

xStart = 1.8
xEnd = 1.9
A_throat = 0.7
nPoints = 100

x = np.linspace(xStart, xEnd, nPoints)

# parabolic function with throat in the middle
Acoeff = -4*(A_throat-1)
Bcoeff = -Acoeff
Ccoeff = 1
z = (x-x[0])/(x[-1]-x[0])
Area = Acoeff*z**2 + Bcoeff*z + Ccoeff

plt.figure()
plt.plot(x, Area)
plt.xlabel('x')
plt.ylabel('A')

with open('nozzle.csv', 'w') as file:
    file.write('x,A\n')
    for i in range(len(x)):
        file.write('%.6f,%.6f\n' %(x[i], Area[i]))

plt.show()



