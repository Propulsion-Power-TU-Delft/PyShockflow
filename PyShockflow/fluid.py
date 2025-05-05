import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


class FluidIdeal():
    """
    Ideal Fluid Class, where thermodynamic properties and transformation are computed with ideal gas laws
    """
    def __init__(self, gmma, Rgas):
        self.gmma = gmma
        self.Rgas = Rgas
    
    def ComputeStaticEnergy_p_rho(self, p, rho):
        return (p / (self.gmma - 1) / rho)
    
    def ComputePressure_rho_e(self, rho, e):
        return (self.gmma-1)*rho*e
    
    def ComputeSoundSpeed_p_rho(self, p, rho):
        return np.sqrt(self.gmma*p/rho)
    
    def ComputeMach_u_p_rho(self, u, p, rho):
        soundSpeed = self.ComputeSoundSpeed_p_rho(p, rho)
        return np.abs(u)/soundSpeed

    def ComputeTemperature_p_rho(self, p, rho):
        return (p/rho)/self.Rgas
    
    def ComputeDensity_p_T(self, p, T):
        return p/self.Rgas/T

    def ComputeEntropy_p_rho(self, p, rho):
        return p/(rho**self.gmma)

    def ComputeFunDerGamma_p_rho(self, p, rho):
        if isinstance(p, np.ndarray): # handle the case when the inputs are arrays
            return 0.5*(self.gmma+1)+np.zeros_like(p)
        else:
            return 0.5*(self.gmma+1)

    def ComputeComprFactorZ_p_rho(self, p, rho):
        if isinstance(p, np.ndarray):
            return 1+np.zeros_like(p)
        else:
            return 1
    
    def ComputeTotalPressure_p_M(self, p, M):
        return p*(1+(self.gmma-1)/2*M**2)**(self.gmma/(self.gmma-1))
    
    def ComputeMach_pt_p(self, pt, p):
        mach = np.sqrt( 2/(self.gmma-1) * ((pt/p)**((self.gmma-1)/self.gmma)-1) )
        return mach
    
    def ComputeTotalTemperature_T_M(self, T, M):
        return T*(1+(self.gmma-1)/2*M**2)
    
    def ComputeTemperature_Tt_M(self, Tt, M):
        return Tt/(1+(self.gmma-1)/2*M**2)
    
    def ComputeInletQuantities(self, pressure, totPressure, totTemperature, direction):
        mach = self.ComputeMach_pt_p(totPressure, pressure)
        temperature = self.ComputeTemperature_Tt_M(totTemperature, mach)
        density = self.ComputeDensity_p_T(pressure, temperature)
        soundSpeed = self.ComputeSoundSpeed_p_rho(pressure, density)
        velocity = mach*soundSpeed*direction
        energy = self.ComputeStaticEnergy_p_rho(pressure, density)
        return density, velocity, energy

    def Compute_gammapv_p_rho(self, p, rho):
        if isinstance(p, np.ndarray):
            gmma_pv = np.zeros_like(p)+self.gmma
        else:
            gmma_pv = self.gmma
        return gmma_pv
    
    
    
    

class FluidReal():
    """
    Real Fluid Class, where thermodynamic properties and transformations are taken from coolprop
    """
    def __init__(self, fluid_name):
        self.fluid_name = fluid_name
    
    def ComputeStaticEnergy_p_rho(self, p, rho):
        e = CP.PropsSI('U', 'P', p, 'D', rho, self.fluid_name)
        return e
    
    def ComputePressure_rho_e(self, rho, e):
        p = CP.PropsSI('P', 'D', rho, 'U', e, self.fluid_name)
        return p
    
    def ComputeSoundSpeed_p_rho(self, p, rho):
        try:
            a = CP.PropsSI("A", "P", p, "D", rho, self.fluid_name)
            return a
        except:
            # two phase region (or close) 
            T = self.ComputeTemperature_p_rho(p, rho)
            try:
                Q = CP.PropsSI("Q", "T", T, "P", p, self.fluid_name)
            except:
                # if the state is very close to saturation line it fails to find the quality -> set artifically to 1
                Q = 1

            # Speed of sound in liquid and vapor phases at the given T and P
            a_liquid = CP.PropsSI("A", "T", T, "Q", 0, self.fluid_name)  # sound speed for liquid phase
            a_vapor = CP.PropsSI("A", "T", T, "Q", 1, self.fluid_name)   # sound speed for vapor phase

            # Calculate weighted speed of sound based on quality
            a = (1 - Q) * a_liquid + Q * a_vapor
            return a
    
    def ComputeMach_u_p_rho(self, u, p, rho):
        soundSpeed = self.ComputeSoundSpeed_p_rho(p, rho)
        return np.abs(u)/soundSpeed
    
    def ComputeTemperature_p_rho(self, p, rho):
        T = CP.PropsSI('T', 'P', p, 'D', rho, self.fluid_name)
        return T
    
    def ComputeDensity_p_T(self, p, T):
        rho = CP.PropsSI('D', 'P', p, 'T', T, self.fluid_name)
        return rho

    def ComputeEntropy_p_rho(self, p, rho):
        s = CP.PropsSI('S', 'P', p, 'D', rho, self.fluid_name)
        return s
    
    def ComputeEntropy_p_T(self, p, T):
        s = CP.PropsSI('S', 'P', p, 'T', T, self.fluid_name)
        return s
    
    def ComputeFunDerGamma_p_rho(self, p, rho):
        try: # if single phase this will work
            G = CP.PropsSI("FUNDAMENTAL_DERIVATIVE_OF_GAS_DYNAMICS", "P", p, "D", rho, self.fluid_name)
            return G
        except: # if close to two phase, we need to do like the speed of sound
            T = self.ComputeTemperature_p_rho(p, rho)
            try:
                Q = CP.PropsSI("Q", "T", T, "P", p, self.fluid_name)
            except:
                # if the state is very close to saturation line it fails to find the quality -> set artifically to 1
                Q = 1

            # G in liquid and vapor phases at the given T
            G_liquid = CP.PropsSI("FUNDAMENTAL_DERIVATIVE_OF_GAS_DYNAMICS", "T", T, "Q", 0, self.fluid_name)  # sound speed for liquid phase
            G_vapor = CP.PropsSI("FUNDAMENTAL_DERIVATIVE_OF_GAS_DYNAMICS", "T", T, "Q", 1, self.fluid_name)   # sound speed for vapor phase

            # Calculate weighted G based on quality
            G = (1 - Q) * G_liquid + Q * G_vapor
            return G
        
    def ComputeComprFactorZ_p_rho(self, p, rho):
        Z = CP.PropsSI('Z', 'P', p, 'D', rho, self.fluid_name)
        return Z

    
    def ComputeInletQuantities(self, pressure, totPressure, totTemperature, direction):
        """The full state must be reconstructed from the quantities given in the arguments.
        The entropy of the static and total state must be the same by definition. This is used to find the temperature.

        Args:
            pressure (float): static pressure
            totPressure (float): total pressure
            totTemperature (float): total temperature
        """
        def compute_function_residual(temperatureGuess):
            entropyStatic = self.ComputeEntropy_p_T(pressure, temperatureGuess)
            entropyTotal = self.ComputeEntropy_p_T(totPressure, totTemperature)
            residual = entropyStatic - entropyTotal
            return residual

        temperature = fsolve(compute_function_residual, totTemperature)[0]
        # if temperature > totTemperature or temperature < totTemperature/2:
        #     temperature = totTemperature*0.9

        
        density = self.ComputeDensity_p_T(pressure, temperature)
        gamma_pv = self.Compute_gammapv_p_rho(pressure, density)
        mach = self.ComputeMach_pt_p_gammapv(totPressure, pressure, gamma_pv)
        soundSpeed = self.ComputeSoundSpeed_p_rho(pressure, density)
        velocity = direction * mach * soundSpeed
        energy = self.ComputeStaticEnergy_p_rho(pressure, density)
        return density, velocity, energy
    
    
    def Compute_gammapv_p_rho(self, p, rho):
        cp = CP.PropsSI("Cpmass", "P", p, "D", rho, self.fluid_name)
        cv = CP.PropsSI("Cvmass", "P", p, "D", rho, self.fluid_name)
        dp_drho_T = CP.PropsSI("d(P)/d(D)|T", "P", p, "D", rho, self.fluid_name)
        dp_dv_T = - rho**2 * dp_drho_T
        gmma_pv = -1/(p*rho) * cp/cv * dp_dv_T
        return gmma_pv
    
    
    def Compute_gammapt_p_T(self, p, T):
        rho = CP.PropsSI("D", "P", p, "T", T, self.fluid_name)
        d_rho_dT_P = CP.PropsSI("d(D)/d(T)|P", "P", p, "T", T, self.fluid_name)
        dv_dT_P = - d_rho_dT_P / (rho**2)
        cp = CP.PropsSI("Cpmass", "P", p, "T", T, self.fluid_name)
        gamma_pT = 1 / (1 - p/cp*dv_dT_P)
        return gamma_pT
    
    
    def ComputeMach_pt_p_gammapv(self, pt, p, gamma_pv):
        """Reference to equation 8.10 Nederstigt MS thesis
        """
        mach = np.sqrt(2/(gamma_pv-1) * ((pt/p)**((gamma_pv-1)/gamma_pv) - 1))
        return mach
        
        
        

            