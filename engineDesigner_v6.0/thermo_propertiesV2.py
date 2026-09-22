import cantera as ct
from CoolProp.CoolProp import PropsSI

class SupercriticalFluid:
    def __init__(self, p_chamber, T_inject, T_initial=None, H_current=None):

        #cantera has the redlich-kwong equations so im using it for the supercritical fluid properties
        self.fluid = ct.Solution('nDodecane_Reitz.yaml', 'nDodecane_RK')
        self.p = p_chamber

        self.fluid.TP = T_inject, self.p
        H_inject = self.fluid.enthalpy_mass

        Tcrit = PropsSI('TCRIT', 'n-Dodecane')
        self.fluid.TP = Tcrit, self.p
        Hcrit = self.fluid.enthalpy_mass

        self.q_barrier = Hcrit - H_inject
        
        #defining thermodynamic state (if we know t_initial, use that, otherwise, use enthalpy)
        if T_initial is not None:
            self.fluid.TP = T_initial, self.p
        else:
            self.fluid.HP = H_current, self.p

        self.film_temp = self.fluid.T
        self.film_H = self.fluid.enthalpy_mass
        self.film_rho = self.fluid.density
        self.film_cp = self.fluid.cp_mass

        #first coolprop fallback in case cantera doesn't have data
        #NIST polynomial fallback for supercritical in case coolprop doesn't have data
        try:
            self.mu = self.fluid.viscosity
            self.k = self.fluid.thermal_conductivity
        except:
            try:
                self.mu = PropsSI('V', 'T', self.film_temp, 'P', self.p, 'n-Dodecane')
                self.k = PropsSI('L', 'T', self.film_temp, 'P', self.p, 'n-Dodecane')
            except:
                self.mu = (1.5e-3) - (1e-6 * self.film_temp) + (1e-9 * self.film_temp**2)
                self.k = 0.14 - (0.0001 * self.film_temp)


        self.film_Pr =  (self.film_cp * self.mu) / self.k

class CoreGas:
    def __init__(self, p_chamber, T_core, U_core, phi):

        #ideal gas assumption
        self.fluid = ct.Solution('gri30.yaml', 'gri30')
        self.fluid.set_equivalence_ratio(phi, 'C3H8', 'O2')
        self.fluid.TP = T_core, p_chamber
        self.fluid.equilibrate('TP')


        self.core_velo = U_core
        self.temp = self.fluid.T
        self.rho = self.fluid.density

        self.cp = self.fluid.cp_mass
        self.cv = self.fluid.cv_mass
        self.gamma = self.cp / self.cv
        self.R = self.cp - self.cv

        self.mu = self.fluid.viscosity
        self.k = self.fluid.thermal_conductivity
            
        self.Pr = (self.cp * self.mu) / self.k



#dumb ahh model
'''
from CoolProp.CoolProp import PropsSI

class film_properties:

    def __init__(self, critical_pressure, film_temp=None, film_enthalpy=None):
        self.fluid = 'dodecane'
        self.critical_pressure = critical_pressure

        #getting critical limits before it goes supercritical
        self.p_crit = PropsSI('PCRIT', self.fluid)
        self.t_crit = PropsSI('TCRIT', self.fluid)


        #determining thermodynamic state
        if film_temp is not None:
            #if temp is known, then calc enthalpy
            self.T_L = film_temp
            self.H_L = PropsSI('H', 'T', self.T_L, 'P', self.p_crit, self.fluid)
        elif film_enthalpy is not None:
            #if enthalpy is known, then calc temp
            self.H_L = film_enthalpy
            self.T_L = PropsSI('T', 'H', self.H_L, 'P', self.p_crit, self.fluid)
        else:
            #error in case neither enthalpy nor temperature are provided
            raise ValueError

        #getting density, viscosity, and specific heat at constant pressure
        #coolprop deals with the supercritical nature of dodecane so no extra work here
        self.rho_L = PropsSI('D', 'H', self.H_L, 'P', self.p_crit, self.fluid)
        self.mu_L = PropsSI('V', 'H', self.H_L, 'P', self.p_crit, self.fluid)
        self.Cp_L = PropsSI('C', 'H', self.H_L, 'P', self.p_crit, self.fluid)


        #finding when either critical limit is passed since surface tension effectively is nonexistent past that
        if self.critical_pressure >= self.p_crit or self.T_L >= self.t_crit:
            #hard code to prevent divison by zero error when above limits
            self.sigma_L = 0.00000000001
        else:
            #if below critical limits, solve for surface tension normally + fallback
            try:
                self.sigma_L = PropsSI('I', 'T', self.T_L, 'P', self.critical_pressure, self.fluid)
            except ValueError:
                #manual fallback value, might change later
                self.sigma_L = 0.00000000001

class gas_properties:
    #core gas properties
    def __init__(self, T_G, critical_pressure, U_G):
        self.fluid = 'dodecane'
        self.T_G = T_G
        self.critical_pressure = critical_pressure
        self.U_G = U_G

        #calculating bulk properties
        self.rho_G = PropsSI('D', 'T', self.T_G, 'P', self.critical_pressure, self.fluid)
        self.mu_G = PropsSI('V', 'T', self.T_G, 'P', self.critical_pressure, self.fluid)
        self.Cp_G = PropsSI('C', 'T', self.T_G, 'P', self.critical_pressure, self.fluid)
        
'''




