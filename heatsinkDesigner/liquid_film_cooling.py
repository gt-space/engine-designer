# assume that coolant remains liquid upon injection

import numpy as np
import scipy.optimize as sci

# Uncomment below and change path if getting "No module named contourDesigner" error
import sys
sys.path.insert(0, '/Users/jacob/Documents/Projects/Research/Film Cooling/engine-designer-master/engineDesigner_v5.0')

from heatsinkDesigner.fuel_props import JetA
from contourDesigner.CEA_properties import ceaToSI
from contourDesigner.CEA_properties import siToCEA

# note : all units except for those used in CEA are standard metric (meters, seconds, Pascals, kg, etc)
class liquid_film_cooling:

    PI = np.pi

    def __init__(self, cea_obj, mdot_gas0, MR, Pc, d_chamber, pressure_chamber, eps, mdot_cool, pressure_orifice_cool, temp_orifice_cool, d_film_orifice, num_film_orifices, orifice_cstar=1):

        self.cea_obj = cea_obj
        self.mdot_gas0 = mdot_gas0 # mass flow rate of non-cooling propellants, kg/s
        self.MR = MR
        self.Pc = Pc # combustion end pressure (Pa)
        self.d_cc = d_chamber
        self.eps = eps # expansion ratio
        self.mdot_cool = mdot_cool
        self.temp_orifice_cool = temp_orifice_cool
        self.pressure0_cool = pressure_orifice_cool
        self.pressure_cc = pressure_chamber # Pascals, this is kinda random for now

        self.area_film_orifice = liquid_film_cooling.PI *(d_film_orifice / 2)**2 * orifice_cstar
        self.rho_inj_cool = JetA.get_rho_l(temp_orifice_cool)
        self.u_inj_cool = self.mdot_cool / (self.rho_inj_cool * self.area_film_orifice * num_film_orifices)
        #print(self.mdot_cool, self.rho_inj_cool, self.area_film_orifice, num_film_orifices, self.u_inj_cool, d_film_orifice, "__________________________________________________U_INJ_COOL TEST__")
        self.mew_inj_l_cool = JetA.get_l_dynamic_viscosity(temp_orifice_cool)

        self.area_cc = np.pi*(self.d_cc/2)**2
        
        self.cp_inj_coolant = JetA.get_c_liquid(self.temp_orifice_cool)
        self.temp_sat_cool = JetA.get_saturation_temp(self.pressure_cc)
        #print(self.temp_sat_cool,"_________________________________________________________________________________________________temp sat cool")
        [ self.molar_mass_gas0, self.cv_gas0, self.cp_gas0, 
          self.mew_gas0, self.prandtal0_gas, self.rho0_gas, self.temp0_chamber,
          self.water_mole_ratio0, self.co2_mole_ratio0 ] = self.get_init_gas_props()
        #print(f"chamber temperature: {self.temp0_chamber}, chamber pressure: {self.pressure_cc}")
        #print (self.molar_mass_gas0, self.cp_gas0, "_____________________________________________________________________________________________cp_gas")

    # https://arc.aiaa.org/doi/pdf/10.2514/6.2004-3360, equation 4
    # simpler model for film cooled length
    # first term finds length before liquid starts boiling, second term finds length after liquid starts boiling
    def solve_fcl_stechman(self):
        sigma = .5 # film stability, this is a minimum value of sigma
        pressure_cc = self.pressure_cc
        mDotL = self.mdot_cool
        temp_recovery = self.temp_orifice_cool
        temp_sat = JetA.get_saturation_temp(pressure_cc)
        # Average coolant cp value between initial and saturated temperature
        cpL = (JetA.get_c_liquid(temp_sat) + self.cp_inj_coolant) / 2
        h_fg = JetA.get_h_fg()
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
    def get_film_cooled_length(self, axial_velocity_fraction):
        PI = liquid_film_cooling.PI
        d_cc = self.d_cc
        area = PI*((d_cc/2)**2)
        circum = PI*d_cc

        pressure = self.pressure_cc
        #print(pressure, "THIS IS THE CHAMBER PRESSUREEEEEEEEEEEEEEEEEEEEEEEEEEEEE")
        temp_sat_cool = JetA.get_saturation_temp(pressure)

        # average of min and max liquid coolant temperatures
        # conservative when considered as a surface temperature
        temp_avg_l_cool = (self.temp_orifice_cool + temp_sat_cool) / 2
        c_l_cool = 2360 #J/kg*K  #JetA.get_c_liquid(temp_avg_l_cool)*3
        #print(c_l_cool, "______________________________________________________________________________________this is the specific heat of the coolant, should be in J/Kg*K")
        rho_l_cool = JetA.get_rho_l(temp_avg_l_cool)
        h_fg = JetA.get_h_fg() # specific heat of vaporization
        mew_l_cool = JetA.get_l_dynamic_viscosity(temp_avg_l_cool) # dynamic viscosity
        #print(mew_l_cool, "_________________________________________________________________________________________this is mew_l_cool")
     
        # these are chamber values without film cooling, calculated by CEA
        pr_gas = self.prandtal0_gas 
        rho_gas = self.rho0_gas
        #print(rho_gas, "____________________________________________________________density of combustion gas")
        cp_gas = self.cp_gas0
        water_mole_ratio = self.water_mole_ratio0
        co2_mole_ratio = self.co2_mole_ratio0
        temp_gas = self.temp0_chamber

        # properties at injection
        u_gas = self.mdot_gas0 / (rho_gas * area)
        u_l_cool = self.u_inj_cool*axial_velocity_fraction #to account for tangential injection

        # density immediately following vaporization
        rho_v_cool = JetA.get_rho_v(temp_sat_cool, pressure)

        # variables to calculate dyamic vsicosity of the coolant vapor
        avg_vapor_mole_flowrate = (self.mdot_cool / 2) / JetA.M # moles/sec of coolant vapor
        avg_cproducts_mole_flowrate = self.mdot_gas0 / self.molar_mass_gas0 # moles/sec of combustion gases
        avg_vapor_mole_ratio = avg_vapor_mole_flowrate / (avg_vapor_mole_flowrate + avg_cproducts_mole_flowrate)
        molecular_g = JetA.M * avg_vapor_mole_ratio + self.molar_mass_gas0 * (1 - avg_vapor_mole_ratio) # average gas molar mass
        mew_v_cool = JetA.get_v_dynamic_viscosity(temp_gas)
        
        # Graham's model
        mew_g = mew_v_cool*avg_vapor_mole_ratio + self.mew_gas0 * (1-avg_vapor_mole_ratio) # weighted average by mole fractions
        #print(rho_v_cool, "rho_v_cool", mew_v_cool, "mew_v_cool", mew_g, "mew_g __________________________________________________ density and mew tests")

        # get h star fg
        hstar_fg = h_fg + (temp_sat_cool-self.temp_orifice_cool)*c_l_cool
        #print(h_fg, hstar_fg, c_l_cool, "_______________________________________________________________________________hstar_fg")

        # Henstock and Hanratty correlation for film thickness (varies depending on liquid film reynolds number)
        interfacial_friction_factor = 0.005 #approx value
        reynold_l_cool = rho_l_cool * u_l_cool * d_cc / mew_l_cool
        interfactial_sheer_stress = 0.5 * interfacial_friction_factor * rho_gas * (u_gas**2)
        friction_velocity = (interfactial_sheer_stress / rho_l_cool)**0.5
        if reynold_l_cool < 1600: #laminar/wavy film
            dimless_thickness = 0.0707 * (reynold_l_cool**0.5)
        else:
            dimless_thickness = 0.051 * (reynold_l_cool**0.875)
        film_thickness = (dimless_thickness * mew_l_cool) / (rho_l_cool * friction_velocity)*100 # the '*100' converts to cm
        #print(film_thickness, "_______________________________________________________________________________________film_thickness")
        # calculating radiative heat trasnfer
        mean_beam_length = film_thickness * 1.8 # formula for radiation exchange between two infinite parallel plates (i.e. the wall and the gas core)
        epsilon = liquid_film_cooling.get_emissitivity(mean_beam_length, temp_gas, pressure*water_mole_ratio, pressure*co2_mole_ratio) #use either 'mean_beam_length' or 'd_cc/4'(original code) for the beam length argument
        boltzman = 5.6703e-8
        rad_heat = boltzman * epsilon * (temp_gas**4-temp_avg_l_cool**4) # radiative heat transfer in J/m^2*s
  
        Kt = 1.6 # assume typical free stream turbulance

        h = self.get_h_g(Kt, hstar_fg, u_gas, u_l_cool, rho_gas, mew_g, molecular_g, pr_gas, cp_gas,temp_gas,temp_sat_cool,rad_heat, d_cc) # h is gas-side heat transfer coefficient
        conv_heat = h * (temp_gas-temp_avg_l_cool) # W/m^2*s
        total_heat_flux = conv_heat + rad_heat
        mdot_v = total_heat_flux / hstar_fg
            
        entrainment_fraction =  self._get_entrainment_fraction(temp_sat_cool, d_cc, rho_l_cool, rho_v_cool, rho_gas, u_l_cool, mew_l_cool, u_gas) / mdot_v

        # coolant flow per circumference available for film cooling
        cool_per_circum = self.mdot_cool*(1-entrainment_fraction)/circum 
        #print(h, "h", epsilon, "epsilon", cool_per_circum, "cool_per_circum", mdot_v, "mdot_v",  entrainment_fraction, "entrainmeent fraction ________________________ fcl tests" )

        # length for which film remains liquid (m)   
        fcl = cool_per_circum/mdot_v 

        return fcl, h, entrainment_fraction

    # used same paper as authors of "A new generalised model for liquid film cooling in rocket combustion chambers":
    # only looks at water vapor and co2 ; apparently these much more significant than other things
    # https://www.sciencedirect.com/science/article/pii/S0010218072800841?pes=vor&utm_source=scopus&getft_integrator=scopus
    @staticmethod
    def get_emissitivity(path_length, temp, pressure_water, pressure_co2):
        lambda_water = np.log(path_length*(pressure_water*(10**-5))) # path length for this is in bar cm
        lambda_co2 = np.log(path_length*(pressure_co2*(10**-5)))
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

    def get_h_g(self, Kt, hstar_fg, u_gas, u_cool, rho_gas, mew_gas, molecular_g, pr_gas,cp_gas,temp_gas,temp_sat_cool, enthalpy_flux_other, diam_cc):
        G_mean = (rho_gas * u_gas) * ((u_gas - u_cool)/u_gas) # free stream gas flux, relative to liquid flow 
        reynoldGas = G_mean * diam_cc / mew_gas
        #print(f"g mean: {G_mean}")

        # # We provide an initial guess for the log of the friction factor
        # log_initial_guess = np.log(0.001)
        # # lambda is Darcy friction factor ; f is fanning friction factor
        # log_f = sci.fsolve(func=liquid_film_cooling.lambda_func_log, x0=log_initial_guess, args=(reynoldGas), xtol=1.49012e-08)[0]
        # lambda_value = np.exp(log_f)
        # # The fanning friction factor (f) is lambda / 4
        # f = lambda_value / 4

        # --- EXPLICIT SOLVER FOR FRICTION FACTOR (Haaland Approximation) ---
        # This replaces the iterative scipy.fsolve to prevent the numerical errors.
        # The Haaland equation is a well-established explicit approximation for the 
        # Colebrook-White equation used to find the friction factor in pipes.
        
        # We must handle the case of a non-positive Reynolds number to prevent log errors.
        if reynoldGas <= 0:
            # This case is not physically realistic for this model but is handled
            # to prevent a crash. A default friction factor is assumed.
            lambda_value = 0.02 
        else:
            # Haaland's explicit equation for smooth pipes (relative roughness = 0)
            term = 6.9 / reynoldGas
            lambda_inv_sqrt = -1.8 * np.log10(term)
            lambda_value = (1 / lambda_inv_sqrt)**2

        # The Fanning friction factor (f) is the Darcy friction factor (lambda) / 4
        f = lambda_value / 4

        #print(f"values of f func (near 0 if numerical solver converged: {liquid_film_cooling.lambda_func(f*4, reynoldGas)}, value of f: {f}")
        stanton_dry = (f/2)/(1.2+11.8*((f/2)**.5)*(pr_gas-1)*(pr_gas**(-1/3)))
        #print(f"stanton dry: {stanton_dry}")
        h_dry = G_mean * cp_gas * stanton_dry * Kt
        
        #get stanton corrected for transpiration
        # --- MODIFIED SOLVER CALL ---
        # Note the new initial guess for log_st
        initial_guess = [np.log(stanton_dry), 0.5] 
        vals, infodict, ier, mesg = sci.fsolve(func=liquid_film_cooling.stanton_F_function_robust, x0=initial_guess, args=(stanton_dry, molecular_g, JetA.M, cp_gas, hstar_fg, temp_gas, temp_sat_cool, enthalpy_flux_other, h_dry), xtol=1.49012e-08, full_output=True)  # Make sure to set this to True

        # --- CONVERT RESULT BACK ---
        if ier == 1: # Check if the solver succeeded
            log_st_solution, F_solution = vals
            st_solution = np.exp(log_st_solution)
            #print(f"Solution found: st = {st_solution}, F = {F_solution}")
        else:
            print(f"Solver did not converge. Message: {mesg}")
            # --- FALLBACK IF SOLVER FAILS ---
            print(f"WARNING: Stanton number solver did not converge. Message: {mesg}")
            print("Using stanton_dry as a fallback value.")
            st_solution = stanton_dry # Use the non-transpiration value as a fallback
      
      
        #print(f"value of st & h func (near 0 if numerical solver converged): {liquid_film_cooling.stanton_F_function(vals, stanton_dry, molecular_g, JetA.M, cp_gas, hstar_fg, temp_gas, temp_sat_cool, enthalpy_flux_other, h_dry)}")
        st_trans = st_solution

        
        #st, F, = liquid_film_cooling.solve_stanton_directly(stanton_dry, molecular_g, JetA.M, cp_gas, hstar_fg, temp_gas, temp_sat_cool, enthalpy_flux_other, h_dry)
        test_h = st_trans * cp_gas * rho_gas * u_gas
        #print(st_trans, test_h, "________________________________________________________________________________________________________This is the test St and F")
        return st_trans * cp_gas * rho_gas * u_gas 

    def _get_entrainment_fraction(self, temp, d_cc, rho_l_cool, rho_v_cool, rho_gas, u_l_cool, mew_l_cool, u_gas):
        # get entrainment fraction
        reynold_l_cool = (rho_l_cool * u_l_cool * d_cc) / mew_l_cool
        #print(rho_l_cool, u_l_cool, d_cc, mew_l_cool, reynold_l_cool, "THIS IS THE TEST VARIABLES---------------")
        deltarho_cool = rho_l_cool - rho_v_cool
        a = 2.31e-4*reynold_l_cool**-0.35
        Em = 1 - (250*np.log(reynold_l_cool)-1265)/reynold_l_cool
        tension = JetA.get_surface_tension() # surface tension
        We = rho_gas*(u_gas**2)*d_cc/(tension*(deltarho_cool/rho_gas)**0.25)
        return Em*np.tanh(a*We**1.25) 

    # solution should make this equal 0
    # def lambda_func(lambda_value, reynold):
    #     # Prevent taking the square root or log of a non-positive number.
    #     # A negative friction factor is not physically possible.
    #     if lambda_value <= 0:
    #         return 1e6 # Return a large number to push the solver away from this invalid region
        
    #     # if lambda_value == 0:
    #     #     lambda_value = -1e20 # avoid dividing by 0

    #     lam = 1.93*np.log(reynold*np.sqrt(lambda_value))-.537-1/np.sqrt(lambda_value)
    #     return lam
    
    # # solution should make this equal 0. Solves for the log of the friction factor to ensure positivity.
    # def lambda_func_log(log_lambda_value, reynold):
    #     lambda_value = np.exp(log_lambda_value)[0]
        
    #     # The Colebrook-White equation
    #     # This calculation is now safe because lambda_value is guaranteed to be positive.
    #     lam = 1.93 * np.log(reynold * np.sqrt(lambda_value)) - 0.537 - (1 / np.sqrt(lambda_value))
    #     return lam
    
    # at solutions, st_eqn, and F_eqn will equal 0
    # equations for including transpiration effects to calculate heat transfer coefficient
    def stanton_F_function(vars, st_dry, molecular_g, molecular_c, cp_gas, hstar_fg, temp_gas, temp_c_sat, enthalpy_flux_other, h_dry):
        [st, F] = vars

        st_eqn = st/st_dry - np.log(1+(F/st)*(molecular_g/molecular_c)**.6)/((F/st)*(molecular_g/molecular_c)**.6)
        F_eqn = F/st - (cp_gas/hstar_fg)*((temp_gas-temp_c_sat)+enthalpy_flux_other/h_dry) #  replacing radiative heat flux term in paper with general enthalpy flux that may include radiative heat flux, wall heat flux, and enthalpy flux of outflowing liquid


        #test calculations:
        C = (cp_gas/hstar_fg)*((temp_gas-temp_c_sat)+enthalpy_flux_other/h_dry)

        return [st_eqn, F_eqn]
    
    # def solve_stanton_directly(st_dry, molecular_g, molecular_c, cp_gas, hstar_fg, temp_gas, temp_c_sat, enthalpy_flux_other, h_dry):
    #     """
    #     Solves for the Stanton number (st) and F directly without a numerical solver.
    #     """
        
    #     # 1. Calculate the constant C, which is equal to the ratio F/st
    #     # This comes from rearranging your second equation: F_eqn.
    #     C = (cp_gas / hstar_fg) * ((temp_gas - temp_c_sat) + enthalpy_flux_other / h_dry)
        
    #     # 2. Define the term inside the logarithm from your first equation (st_eqn)
    #     # Let's call it K.
    #     K = C * (molecular_g / molecular_c)**0.6
        
    #     # 3. Solve for st directly by rearranging your first equation.
    #     # We use np.log1p(K) instead of np.log(1 + K) for better numerical
    #     # stability when K is close to zero. The expression log(1+x)/x approaches 1
    #     # as x approaches 0, and np.log1p helps maintain precision.
    #     if K == 0:
    #         # Handle the edge case to avoid division by zero
    #         log_term = 1.0
    #     else:
    #         log_term = np.log1p(K) / K
            
    #     st = st_dry * log_term
        
    #     # 4. Now that we have st, calculate F using the constant ratio C.
    #     F = C * st
        
    #     return st, F
    
    def stanton_F_function_robust(vars, st_dry, molecular_g, molecular_c, cp_gas, hstar_fg, temp_gas, temp_c_sat, enthalpy_flux_other, h_dry):
        [log_st, F] = vars
        st = np.exp(log_st)  # Transform back to st, guaranteeing it's positive

        # Use a small value to prevent division by zero if st is extremely small
        if st < 1e-12:
            st = 1e-12
            
        # This term appears twice, so calculate it once
        B = (F / st) * (molecular_g / molecular_c)**0.6

        # Guard against the log term becoming invalid
        # If B is close to -1, return a large number to steer the solver away
        if B <= -1.0:
            return [1e6, 1e6] # Return a large error

        # Handle the removable singularity at B=0
        # As B -> 0, log(1+B)/B -> 1
        if abs(B) < 1e-9:
            log_term = 1.0
        else:
            log_term = np.log(1 + B) / B
        
        st_eqn = st / st_dry - log_term
        F_eqn = F / st - (cp_gas / hstar_fg) * ((temp_gas - temp_c_sat) + enthalpy_flux_other / h_dry)

        #  # --- DEBUGGING PRINT 2 ---
        # result = [st_eqn, F_eqn]
        # # Look for a line in your output where one of these types is <class 'numpy.ndarray'>
        # print(f"CALCULATED - Returning: {result} | Types: {type(result[0])}, {type(result[1])}")
        
        return [st_eqn, F_eqn]

    def get_init_gas_props(self):
        # use cea to get properties of combustion products ; these values assume no film cooling
        MR = self.MR
        Pchamber = siToCEA(self.Pc,"pressure")
        mole_fraction_dicts = self.cea_obj.get_SpeciesMoleFractions(
            Pc=Pchamber, MR=MR, eps=self.eps, frozen=0, frozenAtThroat=0, min_fraction=5e-6)
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
            Pc=Pchamber, MR=MR, eps=self.eps, frozen=0)    
        cvGas = ceaToSI(chamberTransport[0], "specific heat")
        cpGas = cvGas + (8.3145 / MGas)  #MGas
        mew_g = ceaToSI(chamberTransport[1], "viscosity") # combustion products viscosity
        prandtl = chamberTransport[3]
        # gas density
        rho_gases = ceaToSI(self.cea_obj.get_Chamber_Density(
            Pc=Pchamber, MR=MR, eps=self.eps), "density")
        temp_chamber = ceaToSI(self.cea_obj.get_Temperatures(
            Pc=self.Pc, MR=self.MR, eps=self.eps, frozen=0, frozenAtThroat=0)[0], "temperature")
        return [ MGas, cvGas, cpGas, mew_g, prandtl, rho_gases, temp_chamber, water_mole_ratio, co2_mole_ratio ]
    
    # self.molar_mass_gas0, self.cv_gas0, self.cp_gas0, 
    #       self.mew_gas0, self.prandtal0_gas, self.rho0_gas, self.temp0_chamber,
    #       self.water_mole_ratio0, self.co2_mole_ratio0 ] = self.get_init_gas_props()