import numpy as np
from utils.fuel_props import JetA
from contourDesigner.CEA_properties import ceaToSI
from contourDesigner.CEA_properties import siToCEA

# assuming that liquid coolant instantly vaporizes upon injection
class gas_film_cooling:
    def __init__(self, beta, S, mdot_cool, mdot_props, num_orifices,cd_orifice,orifice_d,pressure_inj, rho0_comb_gases, pressure_cc, orifice_temp, inj_temp, radii, dz):
        self.beta = beta # coolant injection angle relative to wall, radians
        self.S = S # coolant slot height (m)
        self.mdot_cool = mdot_cool # kg/s
        self.mdot_props = mdot_props # mass flow rate of combustions products, kg/s
        self.pressure_inj = pressure_inj # pressure at the injector, Pa
        self.rho0_comb_gases = rho0_comb_gases # initial density of combustion products, kg/m^3
        self.pressure_cc = pressure_cc
        self.coolant_mass_frac = mdot_cool / mdot_props
        self.orifice_temp = orifice_temp # orifice temperature of coolant (K)
        self.inj_temp = inj_temp # temperature at injection point (K); found by CEA
        self.radii = radii # radii of chamber along different axial locations
        self.dz = dz # distance between each axial location (m)

        self.rho0_cool = JetA.get_rho_v(JetA.get_saturation_temp(self.pressure_inj), self.pressure_inj) # largest possible density value at injector (conservative)
        init_area = np.pi*self.radii[0]**2
        orifice_area = cd_orifice * np.pi*(orifice_d/2)**2
        self.u0_gas = self.mdot_props/ (self.rho0_comb_gases * init_area) # approximate mainstream gas speed as initial speed
        self.u0_cool = self.mdot_cool / (self.rho0_cool * orifice_area)
        self.beta_eff = np.arctan(np.sin(beta)/(np.cos(beta)+self.rho0_comb_gases*self.u0_gas/(self.rho0_cool*self.u0_cool)))
        self.eta= .25 # conservative value for film cooling efficiency

        self.chamber_r = .073

    # estimates wall temperature while coolant is at a different temperature than mainstream gases
    # uses averages / estimations along a specficied target film length
    def get_guess_temp_arr(self, target_fcl, adiabatic_wall_temp):
        end_index = np.floor(target_fcl / self.dz) # index in radii corresponding to target cooled length
        film_cooled_area = self.get_film_cooled_area(end_index)

        avg_temp = (self.orifice_temp+adiabatic_wall_temp)/2 # average chamber temperature (K); probably an overestimate (conservative)
        cp_cool = JetA.get_cp_vapor(avg_temp) # specific heat at constant pressure, J/kgK
        alpha_cool = JetA.get_conductivity(avg_temp)/(JetA.get_rho_v(avg_temp, self.pressure_cc)*cp_cool) # thermal diffusivity
        u_gas = self.u0_gas
        u_cool = self.u0_cool
        avg_L = np.pi*(self.radii[0]+self.radii(end_index))
        h_g = self.get_h_g(u_gas, u_cool, self.eta, cp_cool, alpha_cool, avg_L, target_fcl/2)

        # correlation from H&H
        e_term = np.exp(-h_g/(G_c*cp_cool*self.eta))
        G_c = self.mdot_cool/film_cooled_area # average coolant mass flux (kg/m^2*s)
        temp_term = adiabatic_wall_temp /(adiabatic_wall_temp-self.orifice_temp)
        max_temp = (temp_term-e_term)*(adiabatic_wall_temp-self.orifice_temp) # wall temp at end of film cooled length

        return np.linspace(self.orifice_temp, max_temp, end_index)
    
    def get_target_mdot_cool(self, adiabatic_wall_temp=500,target_fcl=7.5*2.54/100,):
        # end_index = np.floor(target_fcl / self.dz) # index in radii corresponding to target cooled length
        # film_cooled_area = self.get_film_cooled_area(end_index)
        film_cooled_area = (0.001730371398473978/(np.pi*self.radii[1]**2))*(2*np.pi*self.chamber_r)

        avg_temp = (self.orifice_temp+adiabatic_wall_temp)/2 # average chamber temperature (K); probably an overestimate (conservative)
        cp_cool = JetA.get_cp_vapor(avg_temp) # specific heat at constant pressure, J/kgK
        alpha_cool = JetA.get_conductivity(avg_temp)/(JetA.get_rho_v(avg_temp, self.pressure_cc)*cp_cool) # thermal diffusivity
        u_gas = self.u0_gas
        u_cool = self.u0_cool
        # avg_L = np.pi*(self.radii[0]+self.radii[end_index])
        # avg_L = 0.001730371398473978/(np.pi*self.radii[0]**2)
        h_g = self.get_h_g(u_gas, u_cool, cp_cool, alpha_cool, 3, target_fcl/2)

        # correlation from H&H
        # e_term = np.exp(-h_g/(G_c*cp_cool*self.eta))
        temp_term = np.log((adiabatic_wall_temp-400) /(adiabatic_wall_temp-self.orifice_temp))
        G_c = -h_g/(cp_cool*self.eta*temp_term)

    def get_u_inj_cool(self):
        return self.u0_cool

    # https://ntrs.nasa.gov/api/citations/19670008176/downloads/19670008176.pdf, page 8
    def get_h_g(self, u_gas, u_cool, cp_cool, alpha_cool, i, X):
        # convert everything to imperial units for correlation
        u_gas = siToCEA(u_gas, "velocity")
        u_cool = siToCEA(u_cool, "velocity")
        cp_cool = siToCEA(cp_cool, "specific heat")
        alpha_cool = siToCEA(alpha_cool, "conductivity")

        L = 2*np.pi*self.chamber_r/.3048 # circumference at axial location ; meters to feet
        term1 = ((self.S/.3048)*u_gas/alpha_cool)**(1/8) # .3048 is to convert meters to feet
        if u_gas > u_cool:
            term2 = 1+.4*np.arctan(u_gas/u_cool-1)
        else:
            term2 = (u_cool/u_gas)**(1.5*(u_cool/u_gas)-1)
        term3 = np.log(np.cos(.8*self.beta_eff))
        K = .04
        h_g_term = (term3-np.log(self.eta))/(term1*term2)+K
        h_g = h_g_term * 2.205 * self.mdot_props * cp_cool / (L*X/.3048) # 2.205 is to convert kg/s to lbm/s
        return ceaToSI(h_g*.3048**2, "specific heat") # back to standard metric

    # treat each axial segment as cylindrical
    def get_film_cooled_area(self, end_index):
          if end_index <= 0:
              print(f"End index is <= 0 ; try again") # avoid infinite loops
              return 0
          area = 0
          counter = 0
          while True:
            area += np.pi * self.radii[counter]**2 * self.dz
            counter += 1
            if counter == end_index:
                return area