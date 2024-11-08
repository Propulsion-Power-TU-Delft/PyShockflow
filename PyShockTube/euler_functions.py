from NumericalCodes.fluid import FluidIdeal
import numpy as np

def GetPrimitivesFromConservatives(u1, u2, u3, fluid):
        """
        Compute primitive variables from conservative

        Parameters
        -----------

        `u1`: density

        `u2`: density*velocity

        `u3`: density*total energy

        `fluid`: fluid object, ideal or real

        Returns
        -----------

        `rho`: density

        `u`: velocity

        `p`: pressure

        `e`: static energy
        """
        rho = u1
        u = u2/u1
        e = u3/rho - 0.5*u**2
        p = fluid.ComputePressure_rho_e(rho, e)
        return rho, u, p, e


def GetConservativesFromPrimitives(rho, u, p, fluid):
        """
        Compute conservative variables from primitives

        Parameters
        -----------

        `rho`: density

        `u`: velocity

        `p`: pressure

        `fluid`: fluid object, ideal or real

        Returns
        -----------

        `u1`: density

        `u2`: density*velocity

        `u3`: density*total energy

        """
        u1 = rho
        u2 = rho*u
        e = fluid.ComputeStaticEnergy_p_rho(p, rho)
        u3 = rho*(0.5*u**2+e)
        return u1, u2, u3


def EulerFluxFromConservatives(u1, u2, u3, fluid):
        """
        Compute Euler flux vector from conservative variables `u1`, `u2`, `u3`, using a certain `fluid` object. 

        Returns
        --------
        `flux`: flux vector
        """
        rho, u, p, e = GetPrimitivesFromConservatives(u1, u2, u3, fluid)
        et = e+0.5*u**2
        flux = np.zeros(3)
        flux[0] = rho*u
        flux[1] = rho*u**2+p
        flux[2] = u*(rho*et+p)
        return flux