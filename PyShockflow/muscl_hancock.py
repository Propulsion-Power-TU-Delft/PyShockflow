import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
from numpy import sqrt
from PyShockflow.roe_scheme import RoeScheme_Base


class MusclHancock:
    def __init__(self, rhoLL, rhoL, rhoR, rhoRR, uLL, uL, uR, uRR, pLL, pL, pR, pRR, dx):
        """
        Initialize the MUCL-Hancock scheme
        """
        self.gmma = 1.4
        self.rhoLL = rhoLL
        self.rhoL = rhoL
        self.rhoR = rhoR
        self.rhoRR = rhoRR
        self.uLL = uLL
        self.uL = uL
        self.uR = uR
        self.uRR = uRR
        self.pLL = pLL
        self.pL = pL
        self.pR = pR
        self.pRR = pRR
        self.dx = dx
        # self.u1L, self.u2L, self.u3L = self.GetConsFromPrim(rhoL, uL, pL)
        # self.u1R, self.u2R, self.u3R = self.GetConsFromPrim(rhoR, uR, pR)
    
    def ReconstructInterfaceValues(self):
        """
        Reconstruct the interface values using the slopes
        """
        delta_rhoL, delta_rhoR = self.rhoL-self.rhoLL, self.rhoRR-self.rhoR
        delta_uL, delta_uR = self.uL-self.uLL, self.uRR-self.uR
        delta_pL, delta_pR = self.pL-self.pLL, self.pRR-self.pR


        
        self.rhoL, self.rhoR = self.rhoL+delta_rhoL/2, self.rhoR-delta_rhoR/2
        self.uL, self.uR = self.uL+delta_uL/2, self.uR-delta_uR/2
        self.pL, self.pR = self.rhoL+delta_pL/2, self.rhoR-delta_pR/2
    
    def EvolveInterfaceValues(self, dx, dt):
        """
        Evolve the interface values
        """
        # get conservatives
        self.u1L, self.u2L, self.u3L = self.GetConsFromPrim(self.rhoL, self.uL, self.pL)
        self.u1R, self.u2R, self.u3R = self.GetConsFromPrim(self.rhoR, self.uR, self.pR)

        # left and right flux
        fluxL = self.EulerFlux(self.u1L, self.u2L, self.u3L)
        fluxR = self.EulerFlux(self.u1R, self.u2R, self.u3R)

        # evolve the conservatives
        self.u1L = self.u1L+0.5*dt/dx*(fluxL[0]-fluxR[0])
        self.u2L = self.u2L+0.5*dt/dx*(fluxL[1]-fluxR[1])
        self.u3L = self.u3L+0.5*dt/dx*(fluxL[2]-fluxR[2])
        self.u1R = self.u1R+0.5*dt/dx*(fluxL[0]-fluxR[0])
        self.u2R = self.u2R+0.5*dt/dx*(fluxL[1]-fluxR[1])
        self.u3R = self.u3R+0.5*dt/dx*(fluxL[2]-fluxR[2])

        # update the primitives
        self.rhoL, self.uL, self.pL = self.GetPrimitivesFromCons(self.u1L, self.u2L, self.u3L)
        self.rhoR, self.uR, self.pR = self.GetPrimitivesFromCons(self.u1R, self.u2R, self.u3R)
    
    def GetPrimitivesFromCons(self, u1, u2, u3):
        """
        Compute primitive variables from conservative
        """
        rho = u1
        u = u2/u1
        e = u3/rho - 0.5*u**2
        p = (self.gmma-1)*rho*e
        return rho, u, p

    def GetConsFromPrim(self, rho, u, p):
        """
        Compute conservative variables from primitives
        """
        u1 = rho
        u2 = rho*u
        e = p/(self.gmma-1)/rho
        u3 = rho*(0.5*u**2+e)
        return u1, u2, u3
    
        
    def EulerFlux(self, u1, u2, u3):
        """
        Compute Euler equations flux starting from conservative variables. The expressions are given in page 89 Toro book
        """
        flux = np.zeros(3)
        flux[0] = u2
        flux[1] = 0.5*(3-self.gmma)*u2**2/u1+(self.gmma-1)*u3
        flux[2] = self.gmma*u2*u3/u1-0.5*(self.gmma-1)*u2**3/u1**2
        return flux
    
    def ComputeRoeFlux(self):
        """
        Compute the flux at the interface making use of Roe Solver
        """
        roe = RoeScheme_Base(self.rhoL, self.rhoR, self.uL, self.uR, self.pL, self.pR)
        roe.ComputeAveragedVariables()
        roe.ComputeAveragedEigenvalues()
        roe.ComputeAveragedEigenvectors()
        roe.ComputeWaveStrengths()
        flux = roe.ComputeFlux()
        return flux




















