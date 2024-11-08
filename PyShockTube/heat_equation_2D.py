import numpy as np
import matplotlib.pyplot as plt
import math
import os
import pickle


class HeatEquation2D:
    def __init__(self, xgrid, ygrid):
        """
        Initialize the problem with space arrays. Finit difference formulation.
        """
        self.xgrid = xgrid
        self.ygrid = ygrid
        self.nx = xgrid.shape[0]
        self.ny = ygrid.shape[1]
        self.nPoints = self.nx*self.ny
        self.dx = self.xgrid[1, 0]-self.xgrid[0, 0]
        self.dy = self.ygrid[0, 1]-self.ygrid[0, 0]
    
    def SetBoundaryConditions(self, BC_type, *args):
        self.BC_type = BC_type
        if BC_type=='Dirichlet':
            self.Tleft = args[0]
            self.Tbottom = args[1]
            self.Tright = args[2]
            self.Ttop = args[3]
        else:
            raise ValueError('Unknown BC type')
    
    def SetInitialConditions(self, method):
        if method=='Average' and self.BC_type=='Dirichlet':
            Tavg = 0.25*(self.Tleft+self.Tbottom+self.Tright+self.Ttop)
            self.Tgrid = np.zeros_like(self.xgrid)+Tavg
            self.Tvec = self.Tgrid.flatten()
        else:
            raise ValueError('Unknown initialization method or BC type') 
    
    def InitializeLinearSystem(self):
        """
        The system to solve is S*x = f + q, where A comes from the spatial discretization of laplace operator, x is the solution vector,
        f is the source vector (if present), and q is the vector due to the specified boundary conditions on the borders.
        """
        self.S = np.zeros((self.nPoints, self.nPoints))
        self.Q = np.zeros((self.nPoints, 1))
        
        # diagonal terms
        for iPoint in range(0, self.nPoints):
            self.S[iPoint, iPoint] = -2*(self.dx**2+self.dy**2)
        
        for iRow in range(0, self.nPoints):
            i,j = self.GetIJIndex(iRow)
            if i!=0 and i!=self.nx-1 and j!=0 and j!=self.ny-1:
                pl = self.GetPointIndex(i-1, j)  # point index of left neighbor
                pr = self.GetPointIndex(i+1, j)  # point index of right neighbor
                pd = self.GetPointIndex(i, j-1)  # point index of down neighbor
                pu = self.GetPointIndex(i, j+1)  # point index of top neighbor
                self.S[iRow, pr] += self.dy**2
                self.S[iRow, pl] += self.dy**2
                self.S[iRow, pu] += self.dx**2
                self.S[iRow, pd] += self.dx**2
            elif i==0: # we are here on the left border, provide simple equation u = boundary value
                self.S[iRow, iRow] = 1
                self.Q[iRow, 0] = self.Tleft
            elif i==self.nx-1:
                self.S[iRow, iRow] = 1
                self.Q[iRow, 0] = self.Tright
            elif j==0:
                self.S[iRow, iRow] = 1
                self.Q[iRow, 0] = self.Tbottom
            elif j==self.ny-1:
                self.S[iRow, iRow] = 1
                self.Q[iRow, 0] = self.Ttop
        

    def GetIJIndex(self, iPoint):
        """
        Given a certain point, get the order it would have on a 2D array
        """
        i = math.floor(iPoint/self.ny)
        j = iPoint-i*self.ny
        return i,j
    
    def GetPointIndex(self, i, j):
        """
        Given a certain i,j element, get the order it would have on a 1D array
        """
        iPoint = j+i*self.ny
        return iPoint
    
    def SolveSystemDirect(self):
        self.Tvec = np.linalg.inv(self.S)@self.Q
        self.Tgrid = np.reshape(self.Tvec, (self.nx, self.ny))
    

    def PlotSolution(self):
        plt.figure()
        plt.contourf(self.xgrid, self.ygrid, self.Tgrid, levels=30)
        plt.colorbar()
        plt.gca().set_aspect('equal')
    
    def SolveSystemIterative(self, method):
        if method.lower()=='jacobi':
            self.SolveSysJacobi()

    def SolveSysJacobi(self, nIter=100):
        for nIt in range(nIter):
            Res = self.S@self.Tvec-self.Q
            logRes = np.log(np.sum(Res)/self.nPoints)
            for iRow in range(self.nPoints):
                i,j = self.GetIJIndex(iRow)
                if i!=0 and i!=self.nx-1 and j!=0 and j!=self.ny-1:
                    diag = self.S[iRow, iRow]
                    delta = -Res[iRow, 0]/diag
                    self.Tvec[iRow] += delta
            self.Tgrid = np.reshape(self.Tvec, (self.nx, self.ny))
            self.PlotSolution()

        




    

    


        

