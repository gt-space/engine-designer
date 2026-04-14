import cantera as ct

# calculates the heating value of fuel 
def heating_value(T=298, P=2.689e+6):
    gas = ct.Solution('nDodecane_Reitz.yaml', 'nDodecane_RK')
    fuel = 'c12h26' # setting fuel species as nDodecane
    gas.TP = T, P
    gas.set_equivalence_ratio(0.5, 'c12h26:1.0', "O2:1.0") # using gaseous o2, !!may need to account for phase change later
    h1 = gas.enthalpy_mass # enthalpy
    Y_fuel = gas[fuel].Y[0] # mass fraction for stoich to /kg of fuel

    # combustion products
    X_products = {
        'CO2': gas.elemental_mole_fraction('C'),
        'H2O': 0.5 * gas.elemental_mole_fraction('H')
    }

    gas.TPX = None, None, X_products # updating with products
    water = ct.Water() # water calcs for HHV
    water.TQ = 298, 0 # liquid water state (temp, vapor fraction)
    h_liquid = water.h

    water.TQ = 298, 1 # gaseous water state (temp, vapor fraction)
    h_gas = water.h
    Y_H2O = gas["H2O"].Y[0]

    h2 = gas.enthalpy_mass
    LHV = -(h2-h1) / Y_fuel / 1e6 # low heating value
    HHV = -(h2-h1 + (h_liquid - h_gas) * Y_H2O) / Y_fuel / 1e6 # high heating value, accounts for water condensation
                                                               # !!not sure if we need this

    return [LHV, HHV]

# calculates rayleigh line loss based on heating value of fuel and inlet conditions
# inlet mach, heating value, inlet temperature, mdot of air, mdot of fuel
def rayleigh_loss(M, q, T_i, m_air, m_fuel):
    gamma = 1.4 # from seitzman slides, verify
    C_p = 1004 # specific heat of air, J/(kgK)
    q = q * m_air / m_fuel # heating value of fuel
    T_o_i = T_i*(1 + (gamma-1)/2 * M**2) # stagnation temperature, inlet
    T_o_e = T_o_i*(1 + (q*1e6/C_p/T_o_e))

    # TODO: get reference values from table mathmatically
    
    pass



h = heating_value()
print('LHV: ' + str(h[0]) + ' MJ/kg')
print('HHV: ' + str(h[1]) + ' MJ/kg')


# general rayleigh calcs that might be needed
class Rayleigh:
    def __init__(self, roe, p, gamma, M, T):
        self.roe = roe
        self.p = p
        self.gamma = gamma
        self.M1 = M
        self.T = T

    def staticPressure(self):
        # use inital state and M=1 at throat to solve
        
        num = 1 + self.gamma * self.M1**2
        denom = 1 + self.gamma # assume M=1 at exit (choked)
        p_e = num / denom * self.p

        return p_e

    def staticTemperature(self):
        num = (1 + self.gamma * self.M1**2) ** 2
        denom = ((1 + self.gamma) ** 2) * self.M1**2
        T_e = num / denom * self.T

        return T_e

    def stagPressure(self):
        p_static = self.staticPressure()
        num = 1 + (self.gamma-1)/2 * self.M**2
        denom = (self.gamma + 1) / 2
        p_e = p_static * (num/denom)**(self.gamma/(self.gamma-1))
        
        return p_e

    def stagTemperature(self):
        num = (2 * (1 + self.gamma) * self.M**2) * (1 + (self.gamma-1) * self.M**2/2)
        denom = (1 + self.gamma * self.M**2)**2
        T_e = self.T * denom / num

        return T_e

