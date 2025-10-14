from tespy.tools.fluid_properties.wrappers import FluidPropertyWrapper
from .fluid_properties import AbstractState
from .fluid_properties import fluid as fluidclass
from ast import literal_eval
import numpy as np


class FPTespyWrapper(FluidPropertyWrapper):

    def __init__(self, fluid, back_end=None) -> None:
        self.EoS = None
        data = fluid.split('##')

        if len(data) == 2:
            library = data[0]
            fluid_name = data[1]
            self.EoS = AbstractState(fluidclass(library, fluid_name))
        elif len(data) == 3:
            if data[0] == 'pseudofluid':
                library = data[1]
                fluid_name = 'pseudofluid'
                inputs = literal_eval(data[2])
                if library == 'HOGC-PCP-SAFT':
                    n_cmp = inputs['nCmp']
                    cnc = inputs['cnc']
                    groups = np.full((n_cmp, 10), 'None', dtype="U10")
                    groups_occurrences = np.full((n_cmp, 10), -1.0)
                    for i in range(n_cmp):
                        groups[i, 0:len(inputs['Groups'][i])] = inputs['Groups'][i]
                        groups_occurrences[i, 0:len(inputs['Groups Occ'][i])] = inputs['Groups Occ'][i]

                    fluid_object = fluidclass('HOGC-PCP-SAFT', 'pseudofluid', is_pseudofluid=True)
                    fluid_object.set_pseudofluid_groups(n_cmp, cnc, groups, groups_occurrences)
                    fluid_object.all_prop()
                    self.EoS = AbstractState(fluid_object)
            else:
                raise Exception('Fluid name with more than one # divider but it is not a pseudofluid')

        super().__init__(fluid_name, 'FluidProp&' + library)
        self._set_constants()

    def _set_constants(self):
        self._T_min = 200
        self._T_max = 2000
        self._p_min = 0.01
        self._p_max = 100e6
        self._p_crit = self.EoS.p_critical()
        self._T_crit = self.EoS.T_critical()
        self._molar_mass = self.EoS.molar_mass()

    def _is_below_T_critical(self, T):
        return T < self._T_crit

    def _make_p_subcritical(self, p):
        if p > self._p_crit:
            p = self._p_crit * 0.99
        return p

    def get_T_max(self, p):
        return self._T_max

    def isentropic(self, p_1, h_1, p_2):
        return self.h_ps(p_2, self.s_ph(p_1, h_1))

    def T_ph(self, p, h):
        self.EoS.update('Ph', p, h)
        return self.EoS.T()

    def T_ps(self, p, s):
        self.EoS.update('Ps', p, s)
        return self.EoS.T()

    def h_pQ(self, p, Q):
        self.EoS.update('Pq', p, Q)
        return self.EoS.hmass()

    def h_ps(self, p, s):
        self.EoS.update('Ps', p, s)
        return self.EoS.hmass()

    def h_pT(self, p, T):
        self.EoS.update('PT', p, T)
        return self.EoS.hmass()

    def h_QT(self, Q, T):
        self.EoS.update('Tq', T, Q)
        return self.EoS.hmass()

    def s_QT(self, Q, T):
        self.EoS.update('Tq', T, Q)
        return self.EoS.smass()

    def T_sat(self, p):
        p = self._make_p_subcritical(p)
        self.EoS.update('Pq', p, 0)
        return self.EoS.T()

    def p_sat(self, T):
        if T > self._T_crit:
            T = self._T_crit * 0.99

        self.EoS.update('Tq', T, 0)
        return self.EoS.p()

    def Q_ph(self, p, h):
        p = self._make_p_subcritical(p)
        self.EoS.update('Ph', p, h)
        q = self.EoS.Q()
        if 0 < q < 1:
            return q
        T = self.EoS.T()
        self.EoS.update('Pq', p, q)
        if self.EoS.T() == T:
            return q
        else:
            return -1

    def d_ph(self, p, h):
        self.EoS.update('Ph', p, h)
        return self.EoS.rhomass()

    def d_pT(self, p, T):
        self.EoS.update('PT', p, T)
        return self.EoS.rhomass()

    def d_QT(self, Q, T):
        self.EoS.update('Tq', T, Q)
        return self.EoS.rhomass()

    def viscosity_ph(self, p, h):
        self.EoS.update('Ph', p, h)
        return self.EoS.viscosity()

    def viscosity_pT(self, p, T):
        self.EoS.update('PT', p, T)
        return self.EoS.viscosity()

    def s_ph(self, p, h):
        self.EoS.update('Ph', p, h)
        return self.EoS.smass()

    def s_pT(self, p, T):
        self.EoS.update('PT', p, T)
        return self.EoS.smass()
