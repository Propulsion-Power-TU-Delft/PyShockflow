import CoolProp

FPInputs = ['hs', 'Ps', 'PT', 'Pq', 'dh', 'ds','hq','qs']
CPInputs = [CoolProp.HmassSmass_INPUTS, CoolProp.PSmass_INPUTS, CoolProp.PT_INPUTS,
            CoolProp.PQ_INPUTS, CoolProp.DmassHmass_INPUTS, CoolProp.DmassSmass_INPUTS,
            CoolProp.HmassQ_INPUTS,CoolProp.QSmass_INPUTS]

SwapFPInputs = ['Ph', 'Tq', 'Pd', 'Ts']
SwapCPInputs = [CoolProp.HmassP_INPUTS, CoolProp.QT_INPUTS, CoolProp.DmassP_INPUTS, CoolProp.SmassT_INPUTS]

FPProp = ['K', 'ST', 'C', 'CV', 'Mmol']
CPProp = ['L', 'surface_tension', 'CPMASS', 'CVMASS', 'M']

FPLowProp = ['rhomass', 'P', 'T']
CPLowProp = [CoolProp.iDmass, CoolProp.iP, CoolProp.iT]


def translate_coolprop_fluidprop_abstractstate_prop(prop):

    if prop in CPLowProp:
        return FPLowProp[CPLowProp.index(prop)]
    else:
        return prop


def translate_coolprop_fluidprop_abstractstate_input(input_spec, input1, input2):

    if input_spec in CPInputs:
        input_spec = FPInputs[CPInputs.index(input_spec)]
    elif input_spec in SwapCPInputs:
        input_spec = SwapFPInputs[SwapCPInputs.index(input_spec)]
        tmp = input1
        input1 = input2
        input2 = tmp
    else:
        raise Exception('This CoolProp input was not implemented yet. It is very easy to do it, implement it!')

    return input_spec, input1, input2


def translate_fluidprop_coolprop_abstractstate_input(input_spec, input1, input2):

    if input_spec in FPInputs:
        input_spec = CPInputs[FPInputs.index(input_spec)]
    elif input_spec[::-1] in FPInputs:
        input_spec = CPInputs[FPInputs.index(input_spec[::-1])]
        tmp = input1
        input1 = input2
        input2 = tmp
    elif input_spec in SwapFPInputs:
        input_spec = SwapCPInputs[SwapFPInputs.index(input_spec)]
        tmp = input1
        input1 = input2
        input2 = tmp
    else:
        raise Exception('This FluidProp input was not implemented yet. It is very easy to do it, implement it!')

    return input_spec, input1, input2


def translate_fluidprop_coolprop_prop(prop):
    if prop in FPProp:
        return CPProp[FPProp.index(prop)]
    else:
        return prop


def translate_coolprop_fluidprop_prop(prop):
    if prop in CPProp:
        return FPProp[CPProp.index(prop)]
    else:
        return prop


def check_aliases(name):
    if name == 'HFO-1234ze-E':
        return 'R1234ze(E)'

    