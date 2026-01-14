import numpy as np
from scipy.optimize import brentq
from functools import lru_cache

R = 8.31446261815324
class IdealGasFluid:

    """zero of enthalpy and entropy at 298.15 K"""

    def __init__(self, cp_ig_coeff, cp_ig_coeff1000, mmol):

        self.Mmol = mmol
        self.CpIgCoeff = np.asarray(cp_ig_coeff)
        self.CpIgCoeff1000 = np.asarray(cp_ig_coeff1000)
        self.Tvec = np.empty(5)

    def update(self, inputspec, input0, input1):

        self.P = None
        self.T = None

        if inputspec[0] == 'P':
            self.P = input0
        elif inputspec[0] == 'T':
            self.T = input0
        elif inputspec[0] == 'h':
            brentq(self.hloop, 298.15,5000, args=input0)
        elif inputspec[0] == 's':
            brentq(self.sloop, 298.15,5000, args=input0)

        if inputspec[1] == 'P':
            self.P = input1
        elif inputspec[1] == 'T':
            self.T = input1
        elif inputspec[1] == 'h':
            brentq(self.hloop, 298.15,5000, args=input1)
        elif inputspec[0] == 's':
            brentq(self.sloop, 298.15,5000, args=input1)

    def hloop(self, T, h):
        self.T = T
        return  self.hmass() - h

    def sloop(self, T, s):
        self.T = T
        return self.smass() - s

    @lru_cache(maxsize=30)
    def hmass(self):

        """dh = cpdT = C0 + C1T + C2T**2 + C3T**3 + C4T**4 dt"""

        if self.T <= 1000:
            self.Tvec[0] = self.T - 298.15
            self.Tvec[1] = self.T**2 - 298.15**2
            self.Tvec[2] = self.T**3 - 298.15**3
            self.Tvec[3] = self.T**4 - 298.15**4
            self.Tvec[4] = self.T**5 - 298.15**5
            return np.dot(self.Tvec, self.CpIgCoeff)*R/self.Mmol
        else:
            self.Tvec[0] = 1000 - 298.15
            self.Tvec[1] = 1000**2 - 298.15**2
            self.Tvec[2] = 1000**3 - 298.15**3
            self.Tvec[3] = 1000**4 - 298.15**4
            self.Tvec[4] = 1000**5 - 298.15**5
            h1000 = np.dot(self.Tvec, self.CpIgCoeff)
            self.Tvec[0] = self.T - 1000
            self.Tvec[1] = self.T**2 - 1000**2
            self.Tvec[2] = self.T**3 - 1000**3
            self.Tvec[3] = self.T**4 - 1000**4
            self.Tvec[4] = self.T**5 - 1000**5
            return (h1000 + np.dot(self.Tvec, self.CpIgCoeff1000))*R/self.Mmol

    @lru_cache(maxsize=30)
    def smass(self):

        """ds = dq/T = cp/T dT = C0/T + C1 + C2T + C3T**2 + C4T**3 dt """
        if self.T <= 1000:
            self.Tvec[0] = np.log(self.T/298.15)
            self.Tvec[1] = self.T - 298.15
            self.Tvec[2] = self.T**2 - 298.15**2
            self.Tvec[3] = self.T**3 - 298.15**3
            self.Tvec[4] = self.T**4 - 298.15**4
            return np.dot(self.Tvec, self.CpIgCoeff) * R / self.Mmol

        else:
            self.Tvec[0] = np.log(self.T/298.15)
            self.Tvec[1] = 1000 - 298.15
            self.Tvec[2] = 1000**2 - 298.15**2
            self.Tvec[3] = 1000**3 - 298.15**3
            self.Tvec[4] = 1000**4 - 298.15**4
            s1000 = np.dot(self.Tvec, self.CpIgCoeff)
            self.Tvec[0] = np.log(self.T/1000)
            self.Tvec[1] = self.T - 1000
            self.Tvec[2] = self.T**2 - 1000**2
            self.Tvec[3] = self.T**3 - 1000**3
            self.Tvec[4] = self.T**4 - 1000**4
            return (s1000 + np.dot(self.Tvec, self.CpIgCoeff1000))*R/self.Mmol























