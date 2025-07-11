import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
from numpy import sqrt
from PyShockflow.fluid import FluidIdeal
from PyShockflow.euler_functions import *


class RoeScheme_Base:
    def __init__(self, rhoL, rhoR, uL, uR, pL, pR, fluid):
        """
        Roe scheme numerics for ideal gas. Parameters are left and right values of density, velocity and pressure, and the fluid object.
        Formulation based on x-split Riemann Solver in the book by Toro.
        """
        self.rhoL = rhoL
        self.rhoR = rhoR
        self.uL = uL
        self.uR = uR
        self.pL = pL
        self.pR = pR
        self.fluid = fluid
        if isinstance(fluid, FluidIdeal):
            self.gmma = fluid.gmma
        self.eL = fluid.ComputeStaticEnergy_p_rho(pL, rhoL)
        self.eR = fluid.ComputeStaticEnergy_p_rho(pR, rhoR)
        self.htL = self.ComputeTotalEnthalpy(rhoL, uL, pL, self.eL)
        self.htR = self.ComputeTotalEnthalpy(rhoR, uR, pR, self.eR)
        self.u1L, self.u2L, self.u3L = GetConservativesFromPrimitives(rhoL, uL, pL, self.fluid)
        self.u1R, self.u2R, self.u3R = GetConservativesFromPrimitives(rhoR, uR, pR, self.fluid)
        self.aL = self.fluid.ComputeSoundSpeed_p_rho(self.pL, self.rhoL)
        self.aR = self.fluid.ComputeSoundSpeed_p_rho(self.pR, self.rhoR)


    def RoeAVG(self, fL, fR):
        """
        Roe Averaging Operator
        """
        favg = (sqrt(self.rhoL)*fL + sqrt(self.rhoR)*fR)/(sqrt(self.rhoL)+ sqrt(self.rhoR))
        return favg

    
    def ComputeAveragedVariables(self):
        """
        Compute the Roe averaged variables for the 1D Euler equations
        """
        self.rhoAVG = sqrt(self.rhoL*self.rhoR)
        self.uAVG = self.RoeAVG(self.uL, self.uR)
        self.hAVG = self.RoeAVG(self.htL, self.htR)
        self.aAVG = sqrt((self.gmma-1)*(self.hAVG-0.5*self.uAVG**2))
    
    
    def ComputeTotalEnthalpy(self, rho, u, p, e):
        et = 0.5*u**2 + e
        ht = et+p/rho
        return ht
    

    def ComputeAveragedEigenvalues(self):
        """
        Compute eigenvalues of the averaged Jacobian
        """
        self.lambda_vec = np.array([self.uAVG-self.aAVG, 
                                    self.uAVG, 
                                    self.uAVG+self.aAVG])
    

    def ComputeAveragedEigenvectors(self):
        """
        Compute eigenvector matrix of the averaged flux Jacobian
        """
        self.eigenvector_mat = np.zeros((3, 3))
        
        self.eigenvector_mat[0, 0] = 1
        self.eigenvector_mat[1, 0] = self.uAVG-self.aAVG
        self.eigenvector_mat[2, 0] = self.hAVG-self.uAVG*self.aAVG

        self.eigenvector_mat[0, 1] = 1
        self.eigenvector_mat[1, 1] = self.uAVG
        self.eigenvector_mat[2, 1] = 0.5*self.uAVG**2

        self.eigenvector_mat[0, 2] = 1
        self.eigenvector_mat[1, 2] = self.uAVG+self.aAVG
        self.eigenvector_mat[2, 2] = self.hAVG+self.uAVG*self.aAVG
    

    def ComputeWaveStrengths(self):
        """
        Characteristic jumps due to initial conditions
        """
        self.alphas = np.zeros(3)
        self.alphas[0] = 1/2/self.aAVG**2 *(self.pR-self.pL-self.rhoAVG*self.aAVG*(self.uR-self.uL))
        self.alphas[1] = self.rhoR-self.rhoL - (self.pR-self.pL)/self.aAVG**2
        self.alphas[2] = 1/2/self.aAVG**2*(self.pR-self.pL + self.rhoAVG*self.aAVG*(self.uR-self.uL))


    def ComputeLeftRightEigenvalues(self):
        """
        Compute the eigs of left and right values, needed for the entropy fix (Harten-Hyman)
        """
        self.lambda_vecL = np.array([self.uL-self.aL, 
                                    self.uL,
                                    self.uL+self.aL])
        
        self.lambda_vecR = np.array([self.uR-self.aR, 
                                    self.uR,
                                    self.uR+self.aR])
        

    def ComputeFlux(self, entropy_fix=True):
        """
        Compute the Roe flux. The flux is computed for 1D problems.
        """
        fluxL = self.EulerFlux(self.u1L, self.u2L, self.u3L)
        fluxR = self.EulerFlux(self.u1R, self.u2R, self.u3R)
        fluxRoe = 0.5*(fluxL+fluxR)

        # compute the entropy fixed abs eigenvalues
        absEig = np.zeros(3)
        if entropy_fix==False:
            absEig = np.abs(self.lambda_vec)
        else:
            ## Harten-Hymann entropy fix
            self.ComputeLeftRightEigenvalues()
            for k in range(3):
                tmp = np.array([0, self.lambda_vec[k]-self.lambda_vecL[k], self.lambda_vec[k]-self.lambda_vecR[k]])
                delta = np.max(tmp)
                if np.abs(self.lambda_vec[k])<delta:
                    absEig[k] = delta
                else:
                    absEig[k] = np.abs(self.lambda_vec[k])

        for iDim in range(3):
            for jVec in range(3):
                fluxRoe[iDim] -= 0.5*self.alphas[jVec]*absEig[jVec]*self.eigenvector_mat[iDim, jVec]
        
        return fluxRoe
        
    def EulerFlux(self, u1, u2, u3):
        """
        Get the Euler flux starting from conservative variables. 
        """
        flux1D = EulerFluxFromConservatives(u1, u2, u3, self.fluid)
        return flux1D



class RoeScheme_Generalized(RoeScheme_Base):
    """
    Generalised Roe Scheme for real gases, taken from the article 'A simple extension of Roe scheme for real gases', Arabi et al. 
    Journal of Computational Physics 2017. Formulation based on 1D problem.
    """
    def __init__(self, rhoL, rhoR, uL, uR, pL, pR, fluid):
        super().__init__(rhoL, rhoR, uL, uR, pL, pR, fluid)
        self.deltaP = (self.pR-self.pL)
        self.deltaU = (self.uR - self.uL)
        self.deltaRho = (self.rhoR - self.rhoL)
    
    def ComputeAveragedVariables(self):
        """
        Compute the Roe averaged variables for the 1D Euler equations
        """
        self.rhoAVG = sqrt(self.rhoL*self.rhoR)
        self.uAVG = self.RoeAVG(self.uL, self.uR)
        self.hAVG = self.RoeAVG(self.htL, self.htR)
        self.aAVG = self.RoeAVG(self.aL, self.aR)

    def ComputeWaveStrengths(self):
        self.alphas = np.zeros(3)
        self.alphas[0] = 1/2/self.aAVG**2*(self.deltaP+self.rhoAVG*self.aAVG*self.deltaU)
        self.alphas[1] = 1/2/self.aAVG**2*(self.deltaP-self.rhoAVG*self.aAVG*self.deltaU)
        self.alphas[2] = self.deltaRho-self.deltaP/self.aAVG**2
    

    def ComputeAveragedEigenvalues(self):
        self.lambda_vec = np.array([self.uAVG+self.aAVG, 
                                    self.uAVG-self.aAVG,
                                    self.uAVG])
    

    def ComputeLeftRightEigenvalues(self):
        """
        Compute the eigs of left and right values, needed for the entropy fix (Harten-Hyman)
        """
        self.lambda_vecL = np.array([self.uL+self.aL, 
                                    self.uL-self.aL,
                                    self.uL])
        
        self.lambda_vecR = np.array([self.uR+self.aR, 
                                    self.uR-self.aR,
                                    self.uR])
    

    def ComputeFlux(self, entropy_fix=True):
        """
        Assemble the global flux, average + dissipation
        """
        fluxL = self.EulerFlux(self.u1L, self.u2L, self.u3L)
        fluxR = self.EulerFlux(self.u1R, self.u2R, self.u3R)

        self.ComputeDeltaFlux(entropy_fix)
        fluxRoe = 0.5*(fluxL+fluxR) - 0.5*self.deltaF
        return fluxRoe
    

    def ComputeDeltaFlux(self, entropy_fix):
        """
        Compute the correction for the flux. Formulation taken from "A simple extension of Roe's scheme for real gases" by Arabi et al.
        """
        self.deltaF = np.zeros(3)

        # compute the entropy fixed abs eigenvalues
        absEig = np.zeros(3)
        if entropy_fix==False:
            absEig = np.abs(self.lambda_vec)
        else:
            ## Harten-Hymann entropy fix
            self.ComputeLeftRightEigenvalues()
            for k in range(3):
                tmp = np.array([0, self.lambda_vec[k]-self.lambda_vecL[k], self.lambda_vec[k]-self.lambda_vecR[k]])
                delta = np.max(tmp)
                if np.abs(self.lambda_vec[k])<delta:
                    absEig[k] = delta
                else:
                    absEig[k] = np.abs(self.lambda_vec[k])

        self.deltaF[0] = absEig[0]*self.alphas[0] + absEig[1]*self.alphas[1] + absEig[2]*self.alphas[2]
        self.deltaF[1] = (self.uAVG+self.aAVG)*absEig[0]*self.alphas[0] + (self.uAVG-self.aAVG)*absEig[1]*self.alphas[1] + self.uAVG*absEig[2]*self.alphas[2]

        X = (self.rhoR*self.uR*self.htR)-(self.rhoL*self.uL*self.htL)- \
            (self.hAVG+self.uAVG*self.aAVG)*(self.uAVG+self.aAVG)*(1/2/self.aAVG**2*(self.deltaP+self.rhoAVG*self.aAVG*self.deltaU)) - \
            (self.hAVG-self.uAVG*self.aAVG)*(self.uAVG-self.aAVG)*(1/2/self.aAVG**2*(self.deltaP-self.rhoAVG*self.aAVG*self.deltaU))
        
        if (self.uAVG>=0):
            pass
        else:
            X *= -1
        
        self.deltaF[2] = (self.hAVG+self.uAVG*self.aAVG)*absEig[0]*(self.alphas[0]) + \
                         (self.hAVG-self.uAVG*self.aAVG)*absEig[1]*(self.alphas[1]) + X















