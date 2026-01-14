import base_interface
from fluidprop import FPPython as FPPy
import numpy as np
from coolprop import coolprop_functions


class FluidPropFluid(base_interface.Fluid):

    def __init__(self, library, name, arguments):
        super().__init__(library, name)
        self.PrintError = False
        self.IsPseudoFluid = False

        if 'print_error' in arguments.keys():
            self.PrintError = arguments['print_error']
        if 'is_pseudofluid' in arguments.keys():
            self.IsPseudoFluid = arguments['is_pseudofluid']

        self.AreSIunitsSet = False
        if not self.IsPseudoFluid:
            self.get_cmp_cnc()

        FPPy.Init_FluidProp()

    def set_pseudofluid_parameters(self, n_cmp, cnc, sgm_n, sgm_d, epsilon_k,
                                   dipole_moment=0.0, quadrupole_moment=0.0,
                                   is_associating=0.0, association_potential_depth=0.0):

        if self.Library != 'qPCP-SAFT':
            raise Exception('set_pseudofluid_parameters available only with the qPCP-SAFT')
        self.nCmp = n_cmp
        if n_cmp > 1:
            self.cnc = cnc
        else:
            self.cnc = [1.0]

        self.SgmN = sgm_n
        self.SgmD = sgm_d
        self.EpsilonK = epsilon_k

        if dipole_moment == 0.0:
            self.DipoleMoment = [0.0] * self.nCmp
        else:
            self.DipoleMoment = dipole_moment

        if quadrupole_moment == 0.0:
            self.QuadrupoleMoment = [0.0] * self.nCmp
        else:
            self.QuadrupoleMoment = quadrupole_moment

        if is_associating == 0.0:
            self.IsAssociating = [0.0] * self.nCmp
        else:
            self.IsAssociating = is_associating

        if association_potential_depth == 0.0:
            self.AssociationPotentialDepth = [0.0] * self.nCmp
        else:
            self.AssociationPotentialDepth = association_potential_depth

    def set_pseudofluid_groups(self, n_cmp, cnc, group_names, group_occurrences):

        if self.Library != 'HOGC-PCP-SAFT':
            raise Exception('set_pseudofluid_groups available only with the HOGC-PCP-SAFT')

        self.nCmp = n_cmp
        if n_cmp > 1:
            self.cnc = cnc
        else:
            self.cnc = [1.0]

        self.GroupNames = group_names
        self.GroupOccurences = group_occurrences

    def set_fluid(self):
        if self.IsPseudoFluid:
            if self.Library == 'qPCP-SAFT':
                FPPy.SetPseudoFluid(self.Library, self.nCmp, self.cnc, self.SgmN, self.SgmD, self.EpsilonK,
                                    self.DipoleMoment, self.QuadrupoleMoment,
                                    self.IsAssociating, self.AssociationPotentialDepth)
            elif self.Library == 'HOGC-PCP-SAFT':
                FPPy.SetGroups_PseudoFluid(self.nCmp, self.cnc, self.GroupNames, self.GroupOccurences)

        else:
            if self.Library in ['GasMix', 'PCP-SAFT', 'RefProp', 'qPCP-SAFT', 'HOGC-PCP-SAFT']:
                FPPy.SetFluid(self.Library, self.cmp, self.cnc)
            else:
                # joins components with a forward slash between them into a single string
                cmp_str = ['/'.join(self.cmp)]
                FPPy.SetFluid(self.Library, cmp_str, self.cnc)

        if not self.AreSIunitsSet:
            # This setting is very important! FluidProp uses degC, bar and some other non-SI units otherwise!
            FPPy.SetUnits('SI', 'PerMass', 'P,T,h',
                          'Pa,K,J/kg')
            self.AreSIunitsSet = True

    def PropsSI(self, prop, x_str, x, y_str, y):

        self.set_fluid()

        if any([isinstance(x, int), isinstance(x, float), isinstance(y, int), isinstance(y, float)]):
            x = np.array([x])
            y = np.array([y])
        else:
            x = np.array(x)
            y = np.array(y)

        if prop in ['Tcrit', 'Pcrit', 'Tmax', 'Mmol']:
            n = 1
            prop_value = np.zeros(n)

            if prop == 'Tcrit':
                prop_value[0] = FPPy.Tcrit()
            elif prop == 'Pcrit':
                prop_value[0] = FPPy.Pcrit()
            elif prop == 'Tmax':
                prop_value[0] = FPPy.Tmax()
            elif prop == 'Mmol':
                prop_value[0] = FPPy.Mmol()

        else:

            n = len(x)
            prop_value = np.zeros(n)
            # FluidProp only takes pressure and temperature as capital letters.
            # Put all other variable designators to lower case.
            if x_str != 'P' and x_str != 'T':
                x_str = x_str.lower()
            if y_str != 'P' and y_str != 'T':
                y_str = y_str.lower()

            if prop == 'P':
                prop_value = FPPy.Pressure(x_str + y_str, x, y)
            elif prop == 'T':
                prop_value = FPPy.Temperature(x_str + y_str, x, y)
            elif prop == 'D':
                prop_value = FPPy.Density(x_str + y_str, x, y)
            elif prop == 'H':
                prop_value = FPPy.Enthalpy(x_str + y_str, x, y)
            elif prop == 'S':
                prop_value = FPPy.Entropy(x_str + y_str, x, y)
            elif prop == 'Q':
                prop_value = FPPy.VaporQual(x_str + y_str, x, y)
            elif prop == 'U':
                prop_value = FPPy.IntEnergy(x_str + y_str, x, y)
            elif prop == 'C' or prop == 'Cpmass':
                prop_value = FPPy.HeatCapP(x_str + y_str, x, y)
            elif prop == 'A':
                prop_value = FPPy.SoundSpeed(x_str + y_str, x, y)
            elif prop == 'FUNDAMENTAL_DERIVATIVE_OF_GAS_DYNAMICS':
                prop_value = FPPy.Gamma(x_str + y_str, x, y)
            elif prop == 'd(P)/d(D)|T':
                prop_value = FPPy.Chi(x_str + y_str, x, y)
            elif prop == 'd(D)/d(T)|T':
                prop_value = -FPPy.Alpha(x_str + y_str, x, y)*FPPy.Density(x_str + y_str, x, y)
            elif prop == 'Z':
                prop_value = FPPy.Pressure(x_str + y_str, x, y)/(8.3144598/self.Mmol*FPPy.Density(x_str + y_str, x, y))
            elif prop == 'V':
                prop_value = FPPy.Viscosity(x_str + y_str, x, y)
            elif prop == 'K':
                prop_value = FPPy.ThermCond(x_str + y_str, x, y)
            elif prop == 'ST':
                prop_value = FPPy.SurfTens(x_str + y_str, x, y)
            elif prop == 'CV' or prop == 'Cvmass':
                prop_value = FPPy.HeatCapV(x_str + y_str, x, y)
            elif prop == 'Chi':
                prop_value = FPPy.Chi(x_str + y_str, x, y)
            elif prop == 'Kappa':
                prop_value = FPPy.Kappa(x_str + y_str, x, y)

            # The if statement transforms prop_value into a list if it of type np.ndarray
            # This is necessary because some of the used functions of the cycle_package only support certain data types e.g. np.float.

        if len(prop_value) > 1:
            prop_value = list(prop_value)  # Make a list if array np.ndarray dimension larger than one
        else:
            prop_value = prop_value[0]  # Make a np.float64 if only one element

        if not self.PrintError:
            return prop_value

        else:
            if FPPy.GetErrorMsg() != 'No errors':
                print('FluidProp error: ', FPPy.GetErrorMsg(), '\n')
                print('Prop value = -8888.88', '\n')
                prop_value = -8888.88

            return prop_value


class FluidPropAbstractState:

    def __init__(self, fp_fluid, arguments={}):
        self.Fluid = fp_fluid
        self.Fluid.all_prop()
        self.CoolPropLanguage = False
        if 'coolprop_language' in arguments.keys():
            self.CoolPropLanguage = arguments['coolprop_language']

    def update(self, input_spec, input1, input2):
        if self.CoolPropLanguage:
            input_spec, input1,  input2 = coolprop_functions.translate_coolprop_fluidprop_abstractstate_input(input_spec, input1, input2)
        self.Fluid.set_fluid()
        self.AllPropsState = FPPy.AllProps(input_spec, input1, input2)

    def p_critical(self):
        return self.Fluid.Pcrit

    def T_critical(self):
        return self.Fluid.Tcrit

    def molar_mass(self):
        return self.Fluid.Mmol

    def p(self):
        return self.AllPropsState[0]

    def T(self):
        return self.AllPropsState[1]

    def rhomass(self):
        return self.AllPropsState[3]

    def hmass(self):
        return self.AllPropsState[4]

    def smass(self):
        return self.AllPropsState[5]

    def Q(self):
        return self.AllPropsState[7]

    def cvmass(self):
        return self.AllPropsState[10]

    def cpmass(self):
        return self.AllPropsState[11]

    def speed_sound(self):
        return self.AllPropsState[12]

    def fundamental_derivative_of_gas_dynamics(self):
        return self.AllPropsState[22]

    def viscosity(self):
        return self.AllPropsState[23]

    def conductivity(self):
        return self.AllPropsState[24]

    def compressibility_factor(self):
        return self.p() / (8.3144598/self.Fluid.Mmol*self.rhomass()*self.T())

    def first_partial_deriv(self, of, wrt, const):
        if (of == coolprop_functions.CoolProp.iDmass and wrt == coolprop_functions.CoolProp.iP and const == coolprop_functions.CoolProp.iT) or (of == 'rhomass' and wrt == 'P' and const == 'T'):

            ''' here turbosim is requiring (drhodP)_T
            we use chain rule: (drhodP)_T = (drhodv)_T  * (dvdP)_T
            First term (drhodv)_T = d(1/v)dv = -1/v^2 = -rho^2 
            Second term (dvdP)_T  = Kappa/(-1/v) = -Kappa/rho
            Hence:
            (drhodP)_T = Kappa*rho '''

            return self.AllPropsState[21]*self.AllPropsState[3]
