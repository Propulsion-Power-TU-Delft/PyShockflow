R = 8.3144598  # [J/K mol]


def calc_sigma(fluid):
    Tcrit = fluid.PropsSI('Tcrit', [], [], [], [])
    T = 0.7 * Tcrit
    Mmol = fluid.PropsSI('Mmol', [], [], [], [])
    s = fluid.PropsSI('S', 'T', T, 'Q', 1.0)
    sd = fluid.PropsSI('S', 'T', T + 0.01, 'Q', 1.0)

    try:
        sigma = Tcrit / (R / Mmol) * (sd - s) / 0.01
        if sigma == 0:
            sigma = -8888.88

    except Exception as e:
        print('sigma error: ', fluid.library, fluid.name, e)
        sigma = -8888.88

    return sigma


def calc_sigma1(fluid):
    Tcrit = fluid.PropsSI('Tcrit', [], [], [], [])
    Mmol = fluid.PropsSI('Mmol', [], [], [], [])
    cv = fluid.PropsSI('CV', 'P', 0.0001, 'T', Tcrit)
    return 2*cv*Mmol/R
