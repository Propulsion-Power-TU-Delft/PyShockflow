# -*- coding: utf-8 -*-
import numpy as np
from . import fluid_properties as FP
import CoolProp
from scipy.optimize.zeros import brentq


def err_h_sat_L(i, deltaP, deltaH, fluid, Pleft, hleft):
    err_h = hleft + deltaH * i - FP.PropsSI('H', 'P', Pleft + deltaP * i, 'Q', 0.0, fluid)
    return np.round(err_h, 6)


def err_h_sat_V(i, deltaP, deltaH, fluid, Pleft, hleft):
    err_h = hleft + deltaH * i - FP.PropsSI('H', 'P', Pleft + deltaP * i, 'Q', 1.0, fluid)
    return np.round(err_h, 6)


def Set_Prop_Profiles(h_p, P_p, eos, discretization):
    prop_p = PropProfile()
    prop_p.h_p = h_p
    prop_p.P_p = P_p
    prop_p.T_p = np.zeros(discretization)
    prop_p.rho_p = np.zeros(discretization)
    prop_p.q_p = np.zeros(discretization)
    prop_p.mu_p = np.zeros(discretization)
    prop_p.lambda_p = np.zeros(discretization)
    prop_p.Cp_p = np.zeros(discretization)
    prop_p.Data_p = [0] * discretization

    for i in range(discretization):
        eos.update(CoolProp.HmassP_INPUTS, h_p[i], P_p[i])
        prop_p.Data_p[i] = CorrData()
        prop_p.T_p[i] = eos.T()
        prop_p.rho_p[i] = eos.rhomass()
        prop_p.q_p[i] = eos.Q()
        prop_p.mu_p[i] = eos.viscosity()
        prop_p.lambda_p[i] = eos.conductivity()
        prop_p.Cp_p[i] = eos.cpmass()

    return prop_p


def Find_Sat_Conditions(prop_p, eos):
    # check for phase change
    check1 = prop_p.q_p[prop_p.q_p < 1]
    check2 = check1[check1 > 0]

    if check2.size > 0:
        deltaP = (prop_p.P_p[-1] - prop_p.P_p[0])
        deltah = (prop_p.h_p[-1] - prop_p.h_p[0])

        for i in range(0, len(prop_p.h_p) - 1):

            if prop_p.q_p[i] == 0 and prop_p.q_p[i + 1] > 0 and not prop_p.Data_p[i].isSat_L and not prop_p.Data_p[i].isSat_V:
                deltahi_rel = (prop_p.h_p[i] - prop_p.h_p[0]) / deltah
                deltahii_rel = (prop_p.h_p[i + 1] - prop_p.h_p[0]) / deltah
                # found a saturated liquid point

                pos = brentq(err_h_sat_L, deltahi_rel, deltahii_rel,
                             args=(deltaP, deltah, eos.Fluid, prop_p.P_p[0], prop_p.h_p[0]))
                h_satL = prop_p.h_p[0] + (prop_p.h_p[-1] - prop_p.h_p[0]) * pos
                P_satL = prop_p.P_p[0] + (prop_p.P_p[-1] - prop_p.P_p[0]) * pos
                if h_satL == prop_p.h_p[i]:
                    prop_p.Data_p[i].isSat_L = True

                else:
                    eos.update(CoolProp.HmassP_INPUTS, h_satL, P_satL)

                    HTC_data_satL = CorrData()
                    HTC_data_satL.isSat_L = True
                    T_satL = eos.T()
                    rho_satL = eos.rhomass()
                    q_satL = 0.0
                    mu_satL = eos.viscosity()
                    lambda_satL = eos.conductivity()
                    Cp_satL = eos.cpmass()
                    prop_p.h_p = np.insert(prop_p.h_p, i + 1, h_satL)
                    prop_p.P_p = np.insert(prop_p.P_p, i + 1, P_satL)
                    prop_p.T_p = np.insert(prop_p.T_p, i + 1, T_satL)
                    prop_p.rho_p = np.insert(prop_p.rho_p, i + 1, rho_satL)
                    prop_p.q_p = np.insert(prop_p.q_p, i + 1, q_satL)
                    prop_p.mu_p = np.insert(prop_p.mu_p, i + 1, mu_satL)
                    prop_p.lambda_p = np.insert(prop_p.lambda_p, i + 1, lambda_satL)
                    prop_p.Cp_p = np.insert(prop_p.Cp_p, i + 1, Cp_satL)
                    prop_p.Data_p.insert(i + 1, HTC_data_satL)

        # for i in range(0, len(prop_p.h_p) - 1):

            elif 1 > prop_p.q_p[i] >= 0 and (prop_p.q_p[i + 1] == 1 or prop_p.q_p[i + 1] == -1) and \
                    not prop_p.Data_p[i].isSat_L and not prop_p.Data_p[i].isSat_V:
                ### found a saturated vapor point
                deltahi_rel = (prop_p.h_p[i] - prop_p.h_p[0]) / (prop_p.h_p[-1] - prop_p.h_p[0])
                deltahii_rel = (prop_p.h_p[i + 1] - prop_p.h_p[0]) / (prop_p.h_p[-1] - prop_p.h_p[0])
                pos = brentq(err_h_sat_V, deltahi_rel, deltahii_rel,
                             args=(deltaP, deltah, eos.Fluid, prop_p.P_p[0], prop_p.h_p[0]))
                h_satV = prop_p.h_p[0] + deltah * pos
                P_satV = prop_p.P_p[0] + deltaP * pos
                eos.update(CoolProp.HmassP_INPUTS, h_satV, P_satV)

                if h_satV == prop_p.h_p[i]:
                    prop_p.Data_p[i].isSat_V = True

                else:
                    eos.update(CoolProp.HmassP_INPUTS, h_satV, P_satV)
                    HTC_data_satV = CorrData()
                    HTC_data_satV.isSat_V = True

                    T_satV = eos.T()
                    rho_satV = eos.rhomass()
                    q_satV = 1.0
                    mu_satV = eos.viscosity()
                    lambda_satV = eos.conductivity()
                    Cp_satV = eos.cpmass()
                    prop_p.h_p = np.insert(prop_p.h_p, i + 1, h_satV)
                    prop_p.P_p = np.insert(prop_p.P_p, i + 1, P_satV)
                    prop_p.T_p = np.insert(prop_p.T_p, i + 1, T_satV)
                    prop_p.rho_p = np.insert(prop_p.rho_p, i + 1, rho_satV)
                    prop_p.q_p = np.insert(prop_p.q_p, i + 1, q_satV)
                    prop_p.mu_p = np.insert(prop_p.mu_p, i + 1, mu_satV)
                    prop_p.lambda_p = np.insert(prop_p.lambda_p, i + 1, lambda_satV)
                    prop_p.Cp_p = np.insert(prop_p.Cp_p, i + 1, Cp_satV)
                    prop_p.Data_p.insert(i + 1, HTC_data_satV)

    return prop_p


def Mirror_Sat_Conditions(prop_p1, eos1,
                          prop_p2, eos2):
    prop_p1 = Pop_Previous_Mirror_Values(prop_p1)
    prop_p2 = Pop_Previous_Mirror_Values(prop_p2)
    i = 0

    while i < len(prop_p1.h_p) or i < len(prop_p2.h_p):
        if (prop_p1.Data_p[i].isSat_L or prop_p1.Data_p[i].isSat_V) and (
                prop_p2.Data_p[i].isSat_L or prop_p2.Data_p[i].isSat_V):
            deltah_rel1 = (prop_p1.h_p[i] - prop_p1.h_p[0]) / (prop_p1.h_p[-1] - prop_p1.h_p[0])
            deltah_rel2 = (prop_p2.h_p[i] - prop_p2.h_p[0]) / (prop_p2.h_p[-1] - prop_p2.h_p[0])

            if deltah_rel2 <= deltah_rel1:
                prop_p1 = Mirror_Sat_Point(prop_p1, eos1, prop_p2.h_p, i)

            else:
                prop_p2 = Mirror_Sat_Point(prop_p2, eos2, prop_p1.h_p, i)

        elif (prop_p1.Data_p[i].isSat_L or prop_p1.Data_p[i].isSat_V):
            prop_p2 = Mirror_Sat_Point(prop_p2, eos2, prop_p1.h_p, i)

        elif (prop_p2.Data_p[i].isSat_L or prop_p2.Data_p[i].isSat_V):
            prop_p1 = Mirror_Sat_Point(prop_p1, eos1, prop_p2.h_p, i)

        i += 1

    return prop_p1, prop_p2


def Pop_Previous_Mirror_Values(prop_p):
    if len(prop_p.h_p) == len(prop_p.P_p) == len(prop_p.T_p) == len(prop_p.rho_p) == len(prop_p.q_p) == len(
            prop_p.mu_p) == len(prop_p.lambda_p) == len(prop_p.Data_p):

        for i in range(len(prop_p.h_p)):

            if prop_p.Data_p[i].isMirror:
                prop_p.h_p = np.delete(prop_p.h_p, i)
                prop_p.P_p = np.delete(prop_p.h_p, i)
                prop_p.T_p = np.delete(prop_p.h_p, i)
                prop_p.rho_p = np.delete(prop_p.h_p, i)
                prop_p.mu_p = np.delete(prop_p.h_p, i)
                prop_p.lambda_p = np.delete(prop_p.h_p, i)
                prop_p.Cp_p = np.delete(prop_p.Cp_p, i)
                prop_p.Data_p.pop(i)
    else:

        raise Exception('Arrays of different length')

    return prop_p


def Mirror_Sat_Point(prop_p, eos, h_p2, i):
    # profile 2 contains the saturation point, profile 1 needs to find the mirror.

    deltah_rel = (h_p2[i] - h_p2[0]) / (h_p2[-1] - h_p2[0])
    h_Mirror = prop_p.h_p[0] + deltah_rel * (prop_p.h_p[-1] - prop_p.h_p[0])
    P_Mirror = prop_p.P_p[0] + deltah_rel * (prop_p.P_p[-1] - prop_p.P_p[0])
    if prop_p.h_p[i] == h_Mirror:
        prop_p.Data_p[i].isMirror = True

    else:
        eos.update(CoolProp.HmassP_INPUTS, h_Mirror, P_Mirror)

        HTC_data_Mirror = CorrData()
        HTC_data_Mirror.isMirror = True
        T_Mirror = eos.T()
        rho_Mirror = eos.rhomass()
        q_Mirror = eos.Q()
        mu_Mirror = eos.viscosity()
        lambda_Mirror = eos.conductivity()
        Cp_Mirror = eos.cpmass()

        prop_p.h_p = np.insert(prop_p.h_p, i, h_Mirror)
        prop_p.P_p = np.insert(prop_p.P_p, i, P_Mirror)
        prop_p.T_p = np.insert(prop_p.T_p, i, T_Mirror)
        prop_p.rho_p = np.insert(prop_p.rho_p, i, rho_Mirror)
        prop_p.q_p = np.insert(prop_p.q_p, i, q_Mirror)
        prop_p.mu_p = np.insert(prop_p.mu_p, i, mu_Mirror)
        prop_p.lambda_p = np.insert(prop_p.lambda_p, i, lambda_Mirror)
        prop_p.Cp_p = np.insert(prop_p.Cp_p, i, Cp_Mirror)
        prop_p.Data_p.insert(i, HTC_data_Mirror)

    return prop_p


class PropProfile:

    def __init__(self):
        pass


class CorrData:

    def __init__(self):
        self.is_set = False
        self.isSat_L = False
        self.isSat_V = False
        self.isMirror = False
