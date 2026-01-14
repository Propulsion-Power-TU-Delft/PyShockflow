import numpy as np
import fluid_functions


class Fluid:

    def __init__(self, library, name):
        self.Library = library
        self.Name = name
        self.AllPropSet = False

    def get_cmp_cnc(self):

        self.cmp = []  # Mixture components
        self.cnc = []  # Component concentrations
        # ugly sebastian code
        aux_m = self.Name.split('&')  # Separate components list where & occurs
        for aux in aux_m:
            auxx = aux.split('[')  # Separate component name from concentration
            if len(auxx) > 1:  # Executed if more than one component
                self.cmp.append(auxx[0])  # Save mixture component names
                auxxx = float(auxx[1].split(']')[0])
                self.cnc.append(auxxx)  # Save component concentrations
            elif len(auxx) == 1:  # Executed if single component
                self.cmp.append(auxx[0])
                self.cnc.append(1.0)  # Single component concentration is 1
        self.nCmp = len(self.cnc)
        self.cnc = np.array(self.cnc)

    def all_prop(self):

        if not self.AllPropSet:
            self.Tcrit = self.PropsSI('Tcrit', [], [], [], [])
            self.Pcrit = self.PropsSI('Pcrit', [], [], [], [])
            self.Mmol = self.PropsSI('Mmol', [], [], [], [])
            self.sigma = fluid_functions.calc_sigma(self)
            self.sigma1 = fluid_functions.calc_sigma1(self)
            self.AllPropSet = True
