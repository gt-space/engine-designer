# assume that coolant remains liquid upon injection

import numpy as np
import scipy.optimize as sci

# Uncomment below and change path if getting "No module named regenDesigner" error
import sys
sys.path.insert(0, '/Users/atand/OneDrive/Documents/Code stuff/engine-designer/engineDesigner_v5.0')

from regenDesigner.fuel_props import JetA
from contourDesigner.CEA_properties import ceaToSI
from contourDesigner.CEA_properties import siToCEA

# note : all units except for those used in CEA are standard metric (meters, seconds, Pascals, kg, etc)
class liquid_film_cooling:

    PI = np.pi

    def __init__(self, cea_obj, mdot_gas0, MR, Pc, d_chamber, pressure_chamber, eps,
                 mdot_cool, pressure_orifice_cool, temp_orifice_cool, d_film_orifice, num_film_orifices, orifice_cstar, dz):

        self.cea_obj = cea_obj
        self.mdot_gas0 = mdot_gas0 # mass flow rate of non-cooling propellants, kg/s
        self.MR = MR
        self.Pc = Pc # combustion end pressure (Pa)
        self.d_cc = d_chamber
        self.eps = eps # expansion ratio
        self.mdot_cool = mdot_cool
        self.temp_orifice_cool = temp_orifice_cool
        self.pressure0_cool = pressure_orifice_cool
        self.dz = dz
        self.pressure_cc = pressure_chamber # Pascals, this is kinda random for now

        self.area_film_orifice = liquid_film_cooling.PI *(d_film_orifice / 2)**2 * orifice_cstar
        self.rho_inj_cool = JetA.get_rho_l(temp_orifice_cool, self.pressure0_cool)
        self.u_inj_cool = self.mdot_cool / (self.rho_inj_cool * self.area_film_orifice * num_film_orifices)
        self.mew_inj_l_cool = JetA.get_l_dynamic_viscosity(temp_orifice_cool)

        self.d_cc = d_chamber # combustion chamber diameter (m)
        self.area_cc = np.pi*(self.d_cc/2)**2
        
        self.cp_inj_coolant = JetA.get_c_liquid(self.temp_orifice_cool)
        self.temp_sat_cool = JetA.get_saturation_temp(self.pressure_cc)
        [ self.molar_mass_gas0, self.cv_gas0, self.cp_gas0, 
          self.mew_gas0, self.prandtal0_gas, self.rho0_gas, self.temp0_chamber,
          self.water_mole_ratio0, self.co2_mole_ratio0 ] = self.get_init_gas_props()
        print(f"chamber temperature: {self.temp0_chamber}, chamber pressure: {self.pressure_cc}")

    # https://arc.aiaa.org/doi/pdf/10.2514/6.2004-3360, equation 4
    # simpler model for film cooled length
    # first term finds length before liquid starts boiling, second term finds length after liquid starts boiling
    # something might be wrong with this (getting negative h_g values but unsure if this a result of inaccurate inputs)
    def solve_fcl_stechman(self):
        sigma = .5 # film stability, this is a minimum value of sigma
        pressure_cc = self.pressure_cc
        mDotL = self.mdot_cool
        temp_recovery = self.temp_orifice_cool
        temp_sat = JetA.get_saturation_temp(pressure_cc)
        # Average coolant cp value between initial and saturated temperature
        cpL = (JetA.get_c_liquid(temp_sat) + self.cp_inj_coolant) / 2
        h_fg = JetA.get_h_fg(temp_sat)
        cp_gas = ceaToSI(self.cea_obj.get_Chamber_Cp(Pc=siToCEA(self.Pc,"pressure"),MR=self.MR,eps=self.eps),"specific heat")
        hG_no_film = ceaToSI(self.cea_obj.get_Chamber_H(
            Pc=siToCEA(self.Pc, 'pressure'), MR=self.MR, eps=self.eps),'enthalpy')
        lambd = h_fg + (temp_sat-self.temp_orifice_cool)*cpL
        H = cp_gas * (self.temp0_chamber-temp_sat) / lambd
        hG = hG_no_film * (np.log(1+H)/H)
        numerator1 = sigma * mDotL * cpL * (temp_sat - self.u_inj_cool)
        denominator1 = pressure_cc * hG * (temp_recovery - temp_sat)
        numerator2 = sigma * mDotL * h_fg
        denominator2 = pressure_cc * hG * (temp_recovery - temp_sat)
        pre_boil_length = numerator1/denominator1
        return pre_boil_length, pre_boil_length + numerator2 / denominator2

     # https://www.sciencedirect.com/science/article/pii/S0017931012003195#b0085
     # returns total film cooled length for a straight contour
     # analysis approximates values of cp, speed, and temperature to be constant
     # analysis is based on transport of energy for the liquid film
     # heat transfer to the wall is neglected since it is much smaller than heat transfer from mainstream gases (this is conservative --> greater length will be cooled if wall absorbs heat)
    def get_film_cooled_length(self):
        PI = liquid_film_cooling.PI
        d_cc = self.d_cc
        area = PI*((d_cc/2)**2)
        circum = PI*d_cc

        pressure = self.pressure_cc
        temp_sat_cool = JetA.get_saturation_temp(pressure)

        # average of min and max liquid coolant temperatures
        # conservative when considered as a surface temperature
        temp_avg_l_cool = (self.temp_orifice_cool + temp_sat_cool) / 2
        c_l_cool = JetA.get_c_liquid(temp_avg_l_cool)
        rho_l_cool = JetA.get_rho_l(temp_avg_l_cool, pressure)
        h_fg = JetA.get_h_fg(temp_avg_l_cool) # specific heat of vaporization
        mew_l_cool = JetA.get_l_dynamic_viscosity(temp_avg_l_cool) # dynamic viscosity
     
        # these are chamber values without film cooling, calculated by CEA
        pr_gas = self.prandtal0_gas 
        rho_gas = self.rho0_gas
        cp_gas = self.cp_gas0
        water_mole_ratio = self.water_mole_ratio0
        co2_mole_ratio = self.co2_mole_ratio0
        temp_gas = self.temp0_chamber

        # properties at injection
        u_gas = self.mdot_gas0 / (rho_gas * area)
        u_l_cool = self.u_inj_cool

        # density immediately following vaporization
        rho_v_cool = JetA.get_rho_v(temp_sat_cool, pressure)

        avg_vapor_mole_flowrate = (self.mdot_cool / 2) / JetA.M # moles/sec of coolant vapor
        avg_cproducts_mole_flowrate = self.mdot_gas0 / self.molar_mass_gas0 # moles/sec of combustion gases
        avg_vapor_mole_ratio = avg_vapor_mole_flowrate / (avg_vapor_mole_flowrate + avg_cproducts_mole_flowrate)
        molecular_g = JetA.M * avg_vapor_mole_ratio + self.molar_mass_gas0 * (1 - avg_vapor_mole_ratio) # average gas molar mass

        mew_v_cool = JetA.get_v_dynamic_viscosity(temp_gas)
        # Graham's model
        mew_g = mew_v_cool*avg_vapor_mole_ratio + self.mew_gas0 * (1-avg_vapor_mole_ratio) # weighted average by mole fractions

        # get h star fg
        hstar_fg = h_fg + (temp_sat_cool-self.temp_orifice_cool)*c_l_cool

        epsilon = liquid_film_cooling.get_emissitivity(
            d_cc/4, temp_gas, pressure*water_mole_ratio, pressure*co2_mole_ratio)
        boltzman = 5.6703e-8
        rad_heat = boltzman * epsilon * (temp_gas**4-temp_avg_l_cool**4) # radiative heat transfer in J/m^2*s
  
        Kt = 1 # neglect free stream turbulence
        h = self.get_h_g(
            Kt, hstar_fg, u_gas, self.u_inj_cool, rho_gas, mew_g, molecular_g, pr_gas, cp_gas,temp_gas,temp_sat_cool,rad_heat, d_cc) # h is gas-side heat transfer coefficient (J/kgK)
        conv_heat = h * (temp_gas-temp_avg_l_cool) # convective heat flux in J/m^2*s
        total_heat_flux = conv_heat + rad_heat
        mdot_v = total_heat_flux / hstar_fg

        entrainment_fraction =  self._get_entrainment_fraction(
            temp_sat_cool, d_cc, rho_l_cool, rho_v_cool, rho_gas, u_l_cool, mew_l_cool, u_gas) / mdot_v
        # coolant flow per circumference available for film cooling
        cool_per_circum = self.mdot_cool*(1-entrainment_fraction)/circum 
        fcl = cool_per_circum/mdot_v # length for which film remains liquid (m)

        return fcl 

    # used same paper as authors of "A new generalised model for liquid film cooling in rocket combustion chambers":
    # only looks at water vapor and co2 ; apparently these much more significant than other things
    # https://www.sciencedirect.com/science/article/pii/S0010218072800841?pes=vor&utm_source=scopus&getft_integrator=scopus
    @staticmethod
    def get_emissitivity(path_length, temp, pressure_water, pressure_co2):
        lambda_water = np.log(path_length*(pressure_water*10**-5)/100) # path length for this is in bar cm
        lambda_co2 = np.log(path_length*(pressure_co2*10**-5)/100)
        if lambda_water < 0 or lambda_co2 < 0:
            print("Warning: path lengths are too small ; radiative heat being neglected")
            return 0
        tau_water = temp/1000
        tau_co2 = temp/1000

        a0_water = -2.2118 - 1.1987 * tau_water + 0.035596 * tau_water**2
        a1_water = 0.85667 + 0.93048 * tau_water - 0.14391 * tau_water**2
        a2_water = -0.10838 - 0.17156 * tau_water + 0.045915 * tau_water**2

        a0_co2 = -3.3390 + 1.1996 * tau_co2 - 1.0604 * tau_co2**2 + 0.16454 * tau_co2**3
        a1_co2 = 0.90786 + 0.086726 * tau_co2 + 0.13797 * tau_co2**2 - 0.035144 * tau_co2**3
        a2_co2 = -0.15563 -0.10292 * tau_co2 + 0.064443 * tau_co2**2 - 0.014128 * tau_co2**3

        epsilon_water = np.exp(a0_water + a1_water * lambda_water + a2_water * lambda_water ** 2)
        epsilon_co2 = np.exp(a0_co2 + a1_co2 * lambda_co2 + a2_co2 * lambda_co2 ** 2)

        # delta epsilon is a correction factor ; apparently pretty insignificant
        zeta = .5
        avg_lambda = (lambda_water+lambda_co2)/2
        delta_epsilon = (zeta / (10.7 + 101*zeta) - .0089 * zeta**10.4)*avg_lambda**2.76

        return epsilon_water + epsilon_co2 - delta_epsilon

    def get_h_g(self, Kt, hstar_fg, u_gas, u_cool, rho_gas, mew_gas, molecular_g, pr_gas,cp_gas,temp_gas,temp_sat_cool, enthalpy_flux_other, d_cc):
        G_mean = rho_gas * u_gas * (1 - u_cool/u_gas) # free stream gas flux, relative to liquid flow
        reynoldGas = G_mean * d_cc / mew_gas
        print(f"g mean: {G_mean}")

        # lambda is Darcy friction factor ; f is fanning friction factor
        f = sci.fsolve(func=liquid_film_cooling.lambda_func, x0=.001, args=(reynoldGas), xtol=1.49012e-08)[0] / 4
        print(f"values of f func (near 0 if numerical solver converged: {liquid_film_cooling.lambda_func(f*4, reynoldGas)}, value of f: {f}")
        stanton_dry = (f/2)/(1.2+11.8*(f/2)**.5*(pr_gas-1)*pr_gas**(-1/3))
        print(f"stanton dry: {stanton_dry}")
        h_dry = G_mean * cp_gas * stanton_dry * Kt
        
        # get stanton corrected for transpiration
        # st, F, st_dry, molecular_g, molecular_c, cp_gas, hstar_fg, temp_gas, temp_c_sat, qdot_rad, h_dry
        vals = sci.fsolve(func=liquid_film_cooling.stanton_F_function, x0=[stanton_dry, 1], args=(
            stanton_dry, molecular_g, JetA.M, cp_gas, hstar_fg, temp_gas, temp_sat_cool, enthalpy_flux_other, h_dry))
        print(f"value of st & h func (near 0 if numerical solver converged): {liquid_film_cooling.stanton_F_function(vals, stanton_dry, molecular_g, JetA.M, cp_gas, hstar_fg, temp_gas, temp_sat_cool, enthalpy_flux_other, h_dry)}")
        st = vals[0]
        return st * cp_gas * rho_gas * u_gas

    def _get_entrainment_fraction(self, temp, d_cc, rho_l_cool, rho_v_cool, rho_gas, u_l_cool, mew_l_cool, u_gas):
        # get entrainment fraction
        reynold_l_cool = rho_l_cool * u_l_cool * d_cc / mew_l_cool # NO
        deltarho_cool = rho_l_cool - rho_v_cool
        a = 2.31e-4*reynold_l_cool**-0.35
        Em = 1 - (250*np.log(reynold_l_cool)-1265)/reynold_l_cool
        tension = JetA.get_surface_tension(temp) # surface tension
        We = rho_gas*u_gas**2*d_cc/(tension*(deltarho_cool/rho_gas)**0.25)
        return Em*np.tanh(a*We**1.25) 

    # solution should make this equal 0
    def lambda_func(lambda_value, reynold):
        if lambda_value == 0:
            lambda_value = -1e20 # avoid dividing by 0
        return 1.93*np.log(reynold*np.sqrt(lambda_value))-.537-1/np.sqrt(lambda_value)
    
    # at solutions, st_eqn, and F_eqn will equal 0
    def stanton_F_function(vars, st_dry, molecular_g, molecular_c, cp_gas, hstar_fg, temp_gas, temp_c_sat, enthalpy_flux_other, h_dry):
        [st, F] = vars
        st_eqn = st/st_dry - np.log(1+(F/st)*(molecular_g/molecular_c)**.6)/((F/st)*(molecular_g/molecular_c)**.6)
        #  replacing radiative heat flux term in paper with general enthalpy flux that may include radiative heat flux, wall heat flux, and enthalpy flux of outflowing liquid
        F_eqn = F/st - (cp_gas/hstar_fg)*(temp_gas-temp_c_sat+enthalpy_flux_other/h_dry)
        return [st_eqn, F_eqn]

    def get_init_gas_props(self):
        # use cea to get properties of combustion products ; these values assume no film cooling
        MR = self.MR
        Pc = siToCEA(self.Pc,"pressure")
        mole_fraction_dicts = self.cea_obj.get_SpeciesMoleFractions(
            Pc=Pc, MR=MR, eps=self.eps, frozen=0, frozenAtThroat=0, min_fraction=5e-6)
        MGas = 0 # averaged molar mass of combustion products
        for species in mole_fraction_dicts[0].keys():
            MGas += mole_fraction_dicts[1][species][1] * mole_fraction_dicts[0][species]
        MGas /= 1000 # g to kg
        try:
            water_mole_ratio = mole_fraction_dicts[0]['H2O'] / 100
        except KeyError:
            water_mole_ratio = 0
            print("Could not find water mole ratio")
        try:
            co2_mole_ratio = mole_fraction_dicts[0]['*CO2'] / 100
        except KeyError:
            co2_mole_ratio = 0
            print("Could not find CO2 mole ratio")
        chamberTransport = self.cea_obj.get_Chamber_Transport(
            Pc=Pc, MR=MR, eps=self.eps, frozen=0)    
        cvGas = ceaToSI(chamberTransport[0], "specific heat")
        cpGas = cvGas + 8.3145 / MGas
        mew_g = ceaToSI(chamberTransport[1], "viscosity") # combustion products viscosity
        prandtl = chamberTransport[3]
        # gas density
        rho_gases = ceaToSI(self.cea_obj.get_Chamber_Density(
            Pc=Pc, MR=MR, eps=self.eps), "density")
        temp_chamber = ceaToSI(self.cea_obj.get_Temperatures(
            Pc=self.Pc, MR=self.MR, eps=self.eps, frozen=0, frozenAtThroat=0)[0], "temperature")
        return [ MGas, cvGas, cpGas, mew_g, prandtl, rho_gases, temp_chamber, water_mole_ratio, co2_mole_ratio ]