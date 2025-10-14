import base_interface
import numpy as np
from coolprop import coolprop_functions
from lookup_tables.lut import LookUpTable, TableInterpolationError
import copy

InputToLower = ['D', 'Q', 'S']


class LuTFluid(base_interface.Fluid):

    def __init__(self, abstractstate, arguments):

        library = 'LuT'
        super().__init__(library, copy.copy(abstractstate.Fluid.Name))

        self.LuT = []
        self.LuTInputSpecs = []
        self.DefaultLuTIndex = None
        self.add_lut(abstractstate, arguments)

    def add_lut(self, abstractstate, arguments):
        arguments = copy.deepcopy(arguments)
        self.LuT.append(LookUpTable(abstractstate, arguments))
        self.LuTInputSpecs.append(self.LuT[-1].InputSpec)
        if 'default_lut' in arguments.keys():
            if arguments['default_lut']:
                self.DefaultLuTIndex = len(self.LuT) - 1

    def all_prop(self):
        self.Tcrit = self.LuT[0].Tcrit
        self.Pcrit = self.LuT[0].Pcrit
        self.Mmol = self.LuT[0].Mmol
        self.sigma = None
        self.sigma1 = None

    def PropsSI(self, prop, x_str, x, y_str, y):

        prop = coolprop_functions.translate_coolprop_fluidprop_prop(prop)
        if any([isinstance(x, int), isinstance(x, float), isinstance(y, int), isinstance(y, float)]):
            x = np.array([x])
            y = np.array([y])
        else:
            x = np.array(x)
            y = np.array(y)

        if prop in ['Tcrit', 'Pcrit', 'Tmax', 'Mmol']:
            prop_value = np.zeros(1)
            if prop == 'Tcrit':
                prop_value[0] = self.LuT[0].Tcrit
            elif prop == 'Pcrit':
                prop_value[0] = self.LuT[0].Pcrit
            elif prop == 'Mmol':
                prop_value[0] = self.LuT[0].Mmol
        else:        
            n = len(x)
            if x_str in InputToLower:
                x_str = x_str.lower()
            if y_str in InputToLower:
                y_str = y_str.lower()
            prop_value = np.zeros(n)
            inputspec = x_str + y_str
            lutindex = 0
            if inputspec in self.LuTInputSpecs:
                lutindex = self.LuTInputSpecs.index(inputspec)
            
            for i in range(n):
                try:
                    self.LuT[lutindex].update(inputspec, x[i], y[i])
                except TableInterpolationError:
                    if self.DefaultLuTIndex is not None:
                        lutindex = self.DefaultLuTIndex
                        self.LuT[lutindex].update(inputspec, x[i], y[i])
                else:
                    raise

                if prop == 'P':
                    prop_value[i] = self.LuT[lutindex].p()
                elif prop == 'T':
                    prop_value[i] = self.LuT[lutindex].T()
                elif prop == 'D':
                    prop_value[i] = self.LuT[lutindex].rhomass()
                elif prop == 'H':
                    prop_value[i] = self.LuT[lutindex].hmass()
                elif prop == 'S':
                    prop_value[i] = self.LuT[lutindex].smass()
                elif prop == 'Q':
                    prop_value[i] = self.LuT[lutindex].Q()
                elif prop == 'U':
                    prop_value[i] = self.LuT[lutindex].umass()
                elif prop == 'C':
                    prop_value[i] = self.LuT[lutindex].cpmass()
                elif prop == 'A':
                    prop_value[i] = self.LuT[lutindex].sound_speed()
                elif prop == 'V':
                    prop_value[i] = self.LuT[lutindex].viscosity()
                elif prop == 'K':
                    prop_value[i] = self.LuT[lutindex].conductivity()
                elif prop == 'ST':
                    prop_value[i] = self.LuT[lutindex].surface_tension()
                elif prop == 'CV':
                    prop_value[i] = self.LuT[lutindex].cvmass()

        if len(prop_value) > 1:
            prop_value = list(prop_value)  # Make a list if array np.ndarray dimension larger than one
        else:
            prop_value = prop_value[0]  # Make a np.float64 if only one element

        return prop_value


class LuTAbstractState:

    def __init__(self, lut_fluid, arguments={}):
        self.Fluid = lut_fluid
        self.CoolPropLanguage = False
        if 'coolprop_language' in arguments.keys():
            self.CoolPropLanguage = arguments['coolprop_language']

    def update(self, inputspec, input1, input2):
        if self.CoolPropLanguage:
            inputspec, input1,  input2 = (
                coolprop_functions.translate_coolprop_fluidprop_abstractstate_input(inputspec, input1, input2))

        self.ActiveLuTIndex = 0
        if inputspec in self.Fluid.LuTInputSpecs:
            self.ActiveLuTIndex = self.Fluid.LuTInputSpecs.index(inputspec)
        try:
            self.Fluid.LuT[self.ActiveLuTIndex].update(inputspec, input1, input2)
        except TableInterpolationError:
            if self.Fluid.DefaultLuTIndex is not None:
                self.ActiveLuTIndex = self.Fluid.DefaultLuTIndex
                self.Fluid.LuT[self.ActiveLuTIndex].update(inputspec, input1, input2)
            else:
                raise

    def p_critical(self):
        return self.Fluid.LuT[0].p_critical()

    def T_critical(self):
        return self.Fluid.LuT[0].T_critical()

    def molar_mass(self):
        return self.Fluid.LuT[0].molar_mass()

    def p(self):
        return self.Fluid.LuT[self.ActiveLuTIndex].p()

    def T(self):
        return self.Fluid.LuT[self.ActiveLuTIndex].T()

    def rhomass(self):
        return self.Fluid.LuT[self.ActiveLuTIndex].rhomass()

    def hmass(self):
        return self.Fluid.LuT[self.ActiveLuTIndex].hmass()

    def smass(self):
        return self.Fluid.LuT[self.ActiveLuTIndex].smass()

    def Q(self):
        return self.Fluid.LuT[self.ActiveLuTIndex].Q()

    def cvmass(self):
        return self.Fluid.LuT[self.ActiveLuTIndex].cvmass()

    def cpmass(self):
        return self.Fluid.LuT[self.ActiveLuTIndex].cpmass()

    def cp0mass(self):
        return self.Fluid.LuT[self.ActiveLuTIndex].cp0mass()

    def speed_sound(self):
        return self.Fluid.LuT[self.ActiveLuTIndex].speed_sound()

    def fundamental_derivative_of_gas_dynamics(self):
        return self.Fluid.LuT[self.ActiveLuTIndex].fundamental_derivative_of_gas_dynamics()

    def viscosity(self):
        return self.Fluid.LuT[self.ActiveLuTIndex].viscosity()

    def conductivity(self):
        return self.Fluid.LuT[self.ActiveLuTIndex].conductivity()

    def compressibility_factor(self):
        return self.Fluid.LuT[self.ActiveLuTIndex].compressibility_factor()

    def first_partial_deriv(self, of, wrt, const):
        if self.CoolPropLanguage:
            of = coolprop_functions.translate_coolprop_fluidprop_abstractstate_prop(of)
            wrt = coolprop_functions.translate_coolprop_fluidprop_abstractstate_prop(wrt)
            const = coolprop_functions.translate_coolprop_fluidprop_abstractstate_prop(const)
        return self.Fluid.LuT[self.ActiveLuTIndex].first_partial_deriv(of, wrt, const)
