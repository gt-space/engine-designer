from Coolprop.Coolprop import PropsSI

class Liquid_Props:
    def __init__(self, T_L, P_cc):
        self.rho_l = rho_l = PropsSI('D', 'T', T_L, 'P', P_cc, 'dodecane')
        self.mu_l = PropsSI('V', 'T', T_L, 'P', P_cc, 'dodecane')
        self.Cp_l = PropsSI('Cp0mass', 'T', T_L, 'P', P_cc, 'dodecane')
        self.sigma_l = PropsSI('I', 'T', T_L, 'P', P_cc, 'dodecane')
        V = PropsSI('H', 'P', P_cc, 'Q', 1, 'dodecane')
        L = PropsSI('H', 'P', P_cc, 'Q', 0, 'dodecane')
        self.LHV = V-L
        self.T_sat = PropsSI('T', 'P', P_cc, 'Q', 0, 'dodecane')

class Gas_Props:
    def __init__(self, T_L, P_cc):
        self.rho_g = PropsSI('D', 'T', T_L, 'P', P_cc, 'dodecane')
        self.mu_g = PropsSI('V', 'T', T_L, 'P', P_cc, 'dodecane')
        self.Cp_g = PropsSI('Cp0mass', 'T', T_L, 'P', P_cc, 'dodecane')
        self.Pr_g = PropsSI('Prandtl', 'T', T_L, 'P', P_cc, 'dodecane')
