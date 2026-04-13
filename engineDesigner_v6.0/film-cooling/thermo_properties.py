from Coolprop.Coolprop import PropsSI

def get_liquid_props(T_L, P_cc):
    rho_L = PropsSI('D', 'T', T_L, 'P', P_cc, 'dodecane')
    mu_L = PropsSI('V', 'T', T_L, 'P', P_cc, 'dodecane')
    Cp_L = PropsSI('Cp0mass', 'T', T_L, 'P', P_cc, 'dodecane')
    sigma_L = PropsSI('I', 'T', T_L, 'P', P_cc, 'dodecane')
    H_V = PropsSI('H', 'P', P_cc, 'Q', 1, 'dodecane')
    H_L = PropsSI('H', 'P', P_cc, 'Q', 0, 'dodecane')
    LHV = H_V-H_L
    T_sat = PropsSI('T', 'P', P_cc, 'Q', 0, 'dodecane')
    return rho_L, mu_L, Cp_L, sigma_L, LHV, T_sat

def get_gas_props(T_L, P_cc):
    rho_g = PropsSI('D', 'T', T_L, 'P', P_cc, 'dodecane')
    mu_g = PropsSI('V', 'T', T_L, 'P', P_cc, 'dodecane')
    Cp_g = PropsSI('Cp0mass', 'T', T_L, 'P', P_cc, 'dodecane')
    Pr_g = PropsSI('Prandtl', 'T', T_L, 'P', P_cc, 'dodecane')
