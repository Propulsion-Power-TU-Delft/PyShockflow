import matplotlib.pyplot as plt
import numpy as np


class OneDimensional:

    def __init__(self, N):
        self.nPoints = N

    def addCoordinates(self, x):
        self.x = x

    def addDensity(self, rho):
        self.rho = rho

    def addBlockage(self, b):
        self.b = b

    def addArea(self, A):
        self.A = A

    def plotData(self, folder=None):
        plt.figure()
        plt.plot(self.x, self.rho, label='Density')
        plt.plot(self.x, self.b, label='Blockage')
        # plt.plot(self.x, self.A, label='Area')
        plt.xlabel('x')
        plt.legend()
        plt.grid(alpha=0.3)
        if folder is not None:
            plt.savefig(folder + '/data.pdf', bbox_inches='tight')

    def setBC(self, value_0, value_1):
        self.u = np.zeros_like(self.x)
        self.u[0] = value_0
        self.u[-1] = value_1

    def solveContinuity(self, method):
        u = np.zeros_like(self.u)
        u[0] = self.u[0]
        u[-1] = self.u[-1]

        if method == 'original':
            for i in range(1, self.nPoints-1):
                u[i] = (self.rho[i-1]*u[i-1]*self.b[i-1])/(self.b[i]*self.rho[i])
        elif method == 'separated':
            dbdx = np.gradient(self.b, self.x)

            plt.figure()
            plt.plot(self.x, self.b, label='b')
            plt.plot(self.x, dbdx, label='db/dx')
            plt.grid(alpha=0.3)
            plt.xlabel('x')
            plt.legend()
            # plt.savefig('pictures/blockage.pdf', bbox_inches='tight')
            for i in range(1, self.nPoints-1):
                dx = self.x[i]-self.x[i-1]
                u[i] = (self.rho[i-1]*u[i-1]) / (self.rho[i]*(1+dbdx[i]/self.b[i]*dx))
        elif method == 'FVM':
            for i in range(1, self.nPoints-1):
                u[i] = self.rho[i-1]*u[i-1]*(self.b[i-1]+self.b[i])/(self.b[i]+self.b[i+1])
        else:
            raise ValueError('Unknown method!')




        return u

