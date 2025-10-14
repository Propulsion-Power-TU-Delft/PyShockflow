from CoolProp.HumidAirProp import HAPropsSI
import CoolProp


def get_absolute_humidity(T,P,R):

    return HAPropsSI('W', 'T', T, 'P', P, 'R', R)


def translate_prop(prop):
    if prop == 'V':
        return 'mu'
    if prop == 'D':
        return 'Vha'
    if prop == 'H':
        return 'Hha'
    if prop == 'S':
        return 'Sha'
    if prop == 'C':
        return 'Cha'
    if prop == 'CV':
        return 'CVha'
    else:
        return prop


def translate_coolprop_abstractstate_inputspec(inputspec):

    if inputspec == CoolProp.HmassSmass_INPUTS:
        inputspec1 = 'Hha'
        inputspec2 = 'Sha'

    elif inputspec == CoolProp.HmassP_INPUTS:
        inputspec1 = 'Hha'
        inputspec2 = 'P'

    elif inputspec == CoolProp.PSmass_INPUTS:
        inputspec1 = 'P'
        inputspec2 = 'Sha'

    elif inputspec == CoolProp.PT_INPUTS:
        inputspec1 = 'P'
        inputspec2 = 'T'

    elif inputspec == CoolProp.QT_INPUTS:
        inputspec1 = 'Q'
        inputspec2 = 'T'

    elif inputspec == CoolProp.PQ_INPUTS:
        inputspec1 = 'P'
        inputspec2 = 'Q'

    elif inputspec == CoolProp.DmassHmass_INPUTS:
        inputspec1 = 'D'
        inputspec2 = 'Hha'

    elif inputspec == CoolProp.DmassSmass_INPUTS:
        inputspec1 = 'D'
        inputspec2 = 'Sha'

    else:
        raise Exception('This CoolProp input was not implemented yet. It is very easy to do it, implement it!')

    return inputspec1, inputspec2

