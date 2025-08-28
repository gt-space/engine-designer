import cantera as ct
import numpy as np
import scipy.optimize as sci

class CanteraProperties:

    def __init__(self,P_inj,P_amb,con_rat, LOX_inj_temp, JetA_inj_temp,gas_init_temp,MR_weight,max_iters_throat=1_000,max_iters_per_loc=100):

        self.P_inj = P_inj
        self.P_amb = P_amb
        self.con_rat = con_rat
        self.MR = MR_weight * 166 / 32 # ratio of moles
        self.max_iters_throat = max_iters_throat
        self.max_iters_per_loc = max_iters_per_loc

        self.contour_defined = False

        # The below fields are initialized when set_contour is called
        self.Mach_arr = None
        self.area_arr = None
        self.throat_area = None
        self.pressure_arr = None
        self.temp_arr = None
        self.rho_arr = None
        self.ivac_arr = None
        self.cf_arr = None
        self.isp_arr = None
        self.rho_arr = None
        self.h_arr = None
        self.e_arr = None
        self.molar_weight_arr = None
        self.cp_arr = None
        self.gamma_arr = None
        self.sound_speed = None
        self.viscosity = None
        self.thermal_conductivity = None
        self.prandtl_arr = None
        self.velocity_arr = None

        # extract all species in the NASA database
        full_species = {S.name: S for S in ct.Species.list_from_file('nasa_gas.yaml')}

        LOX = ct.Solution('Cantera YAML Files/LOX.yaml','liquid_oxygen')
        JetA = ct.Solution('Cantera YAML Files/JetA.yaml','jet_a')

        # LOX = ct.Solution(thermo='ideal-gas',species=[full_species['O2']])
        # JetA = ct.Solution(thermo='ideal-gas',species=[full_species['Jet-A(g)']])

        # extract only the relevant species
        species = [full_species[S] for S in (
            'CH4', 'CO', 'CO2', 'C2H4', 'H2', 'H2O', 'C'
            )]
        
       
        # only the equilibrium species
        # species = [full_species[S] for S in (
        #     'CH4', 'CO', 'CO2', 'C2H4', 'H2', 'H2O'
        # )]
        self.gas = ct.Solution(thermo='ideal-gas', species=species)
        # spec = ct.Species.list_from_file("Cantera YAML Files/CombProducts.yaml")
        # self.gas = ct.Solution(thermo='ideal-gas', kinetics='gas',species=spec)

        # self.gas = ct.Solution('Cantera YAML Files/CombProducts.yaml')
        self.mixture = ct.Mixture([(LOX,self.MR/(1+self.MR)),(JetA,1/(1+self.MR)),(self.gas,1e-5)]) # does quantity matter?
        
        # https://rocketcea.readthedocs.io/en/latest/finite_area_comb.html
        self.Pc = P_inj/(1+.54/(con_rat**2.2))

        chamber_pressure = 1797626.2626262628 # TEST

        # Find throat properties
        self.h_chamber, self.chamber_gamma, self.chamber_entropy = self.get_chamber_props(chamber_pressure, LOX, JetA, LOX_inj_temp,JetA_inj_temp,P_inj,P_inj,gas_init_temp)
        self.throat_pressure, self.throat_MWt, self.throat_gamma = self.get_throat_props(chamber_pressure, self.h_chamber)

        # set expansion ratio
        k = self.throat_gamma # assume "frozen" at throat for expansion ratio calculation (reaction stops at the throat)
        self.exit_Mach = np.sqrt((2/(k-1))*(Pc/P_amb)^((k-1)/k)-1)
        self.exp_rat  (1/self.exit_Mach)*((1+(k+1)/2*self.exit_Mach^2)/((k+1)/2))^((k+1)/(2*(k-1)))

    def get_Throat_PcOvPe(self):
        return self.throat_pressure / self.P_amb

    def get_eps(self):
        return self.exp_rat

    def get_exit_Mach(self):
        return self.exit_Mach

    def get_throat_MWt(self):
        return self.throat_MWt

    def get_throat_gamma(self):
        return self.throat_gamma

    def get_contour_defined(self):
        return self.contour_defined
    
    '''
    The following methods return a numpy array of engine properties at each axial location
    Once the throat radii are determined, set_contour should be called. Then, the functions for getting properties may be called.
    If a property getter function is called before set_contour is called, an exception will be raised.
    '''

    # Takes a numpy array of radii at each axial location
    def set_contour(self,contour):
        self.area_arr = .5 * np.pi * contour**2
        self.throat_area = np.min(self.area_arr)
        self.set_props()
        self.contour_defined = True

    # Sets properties throughout engine (called internally from set_contour)
    def set_props(self):
        self.pressure_arr = np.array()
        self.temp_arr = np.array()
        self.rho_arr = np.array()
        self.ivac_arr = np.array()
        self.cf_arr = np.array()
        self.isp_arr = np.array()
        self.h_arr = np.array()
        self.e_arr = np.array()
        self.molar_weight_arr = np.array()
        self.cp_arr = np.array()
        self.gamma_arr = np.array()
        self.sound_speed_arr = np.array()
        self.viscosity_arr = np.array()
        self.thermal_conductivity_arr = np.array()
        self.prandtl_arr = np.array()

        '''
        Functions for determining ivac (specific impulse in vaccum), cf (thrust coefficient), isp (specific impulse),
        gamma, viscosity, thermal conductivity, and prandtal number
        '''

        def compute_ivac(Isp,temp,R,Mwt):
            return Isp + temp * R / (Isp * Mwt)
        
        def compute_cstar(gamma, temperature, molecular_weight):
            return (
                np.sqrt(ct.gas_constant * temperature / (molecular_weight * gamma)) *
                np.power(2 / (gamma + 1), -(gamma + 1) / (2*(gamma - 1)))
                )

        def compute_isp(enthalpy):
            return np.sqrt(2 * (self.h_chamber - enthalpy))

        # def compute_viscosity():

        # def compute_thermal_conductivity():

        # def compute_prandtl():

        for i in range(np.size(self.area_arr)):
            mixture = self.mixture
            gas = self.gas
            pressure, gamma, cp = self.get_SP_gamma_cp(self.area_arr[i])

            # Solve for equilibrium at axial location
            gas.SP = self.chamber_entropy, pressure
            mixture.equilibrate('SP')
            gas()

            # Extract properties at axial location
            Isp = compute_isp(gas.enthalpy_mass)
            self.pressure_arr.append(pressure)
            self.temp_arr.append(gas.T)
            self.rho_arr.append(gas.density_mass)
            self.ivac_arr.append(compute_ivac(Isp,gas.T,gas.gas_constant,gas.mean_molecular_weight))
            self.cf_arr.append(Isp / compute_cstar(gamma,gas.T,gas.mean_molecular_weight/1000))
            self.isp_arr.append(Isp)
            self.h_arr.append(gas.enthalpy_mass)
            self.e_arr.append(gas.int_energy_mass)
            self.molar_weight_arr.append(gas.mean_molecular_weight) / 1000 # kg/kmol to kg/mol
            self.cp_arr.append(cp)
            self.gamma_arr.append(gamma)
            self.sound_speed_arr.append(gas.sound_speed)
            self.viscosity_arr.append(0) # TBD
            self.thermal_conductivity_arr.append(0) # TBD
            self.prandtal_arr.append(0) # TBD

        self.velocity_arr = self.Mach_arr * self.sound_speed_arr
    
    def get_SP_gamma_cp(self,area):
        Ae_ov_At = area / self.throat_area

        def iterate(pressure,pinf_pe):
            self.gas.SP = self.chamber_entropy,pressure
            self.mixture.equilibrate('SP')
            self.gas()

            gamma, cp = CanteraProperties._compute_gamma_cp(self.gas)

            speed = np.sqrt(2 * (self.h_chamber - cp))
            speed_sound = self.gas.sound_speed

            Ae_ov_At_est = self.gas.T / (self.gas.P * speed * self.gas.mean_molecular_weight) / self.mdot_tot
            dlogp_dlogA = gamma * speed**2 / (speed**2 - speed_sound**2)
            residual = dlogp_dlogA * (np.log(Ae_ov_At) - np.log(Ae_ov_At_est))
            log_pinf_pe = np.log(pinf_pe/pressure) + residual

            next_pinf_pe = np.exp(log_pinf_pe)
            next_pressure = pinf_pe / next_pinf_pe
    
            return abs(residual) < 4e-5, next_pressure, next_pinf_pe, gamma, cp
        
        converged = False
        pinf_pe = np.exp(self.gamma_throat + 1.4 * np.log(Ae_ov_At))
        pressure = self.chamber_pressure / pinf_pe
        i = 0
        while not converged and i < self.max_iters_per_loc:
            converged, pressure, p_inf_pe, gamma, cp = iterate(pressure,p_inf_pe)
        
        if not converged:
            print(f'Warning: solver did not converge for area={area}')
        
        return pressure, gamma, cp
    
    def get_pressure_arr(self):
        if not self.area_arr:
            raise Exception('Contour must be initialized before property pressure is determined')
        return self.pressure_arr

    def get_temperatures(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property temperature is determined')
        return self.temp_arr
        
    def get_PinjOvP(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property PinjOvP is determined')
        return self.P_inj / self.pressure_arr

    def get_AeOvAt(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property AeOvAt is determined')
        return self.area_arr / self.throat_area
        
    def get_CF(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property CF is determined')
        return self.cf_arr
                

    def get_IVAC(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property IVAC is determined')
        return self.ivac_arr
        

    def get_ISp(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property ISp is determined')
        return self.isp_arr
        
    def get_rho(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property density is determined')
        return self.rho_arr        

    def get_h(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property enthalpy is determined')
        return self.h_arr

    def get_e(self): # return vector (internal energy)
        if not self.area_arr:
            raise Exception('Contour must be initialized before property internal energy is determined')
        return self.e_arr
        
    def get_MWt(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property molar weight is determined')
        return self.molar_weight_arr

    def get_Cp(self): # return vector (specific heat capacity)
        if not self.area_arr:
            raise Exception('Contour must be initialized before property cp is determined')
        return self.cp_arr

    def get_gammas(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property gamma is determined')
        return self.gamma_arr
    
    def get_sound_speeds(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property sound speed is determined')
        return self.sound_speed_arr

    def get_viscosity(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property viscosity is determined')
        return self.viscosity_arr

    def get_conductivity(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property thermal conductivity is determined')
        return self.thermal_conductivity_arr

    def get_prandtl(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property prandtl number is determined')
        return self.prandtl_arr

    def get_velocity_contour(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property velocity is determined')
        return self.velocity_arr
    
    def get_chamber_props(self,chamber_pressure,LOX,JetA,LOX_inj_temp, JetA_inj_temp,LOX_inj_p,JetA_inj_p,gas_init_temp):
        # Find solutions with initial enthalpy of each component

        '''
        CEA Output:
        REACTANT                    WT FRACTION                  ENERGY      TEMP
                                    (SEE note)                  KJ/KG-MOL      K
        FUEL        Jet-A(L)                     1.0000000   -303200.254    298.150
        OXIDANT     O2(L)                        1.0000000    -12978.768     90.180
        '''

        chamber_pressure=self.Pc

        print(f'LOX: {LOX.HP}, JetA: {JetA.HP}, gas: {self.gas.HP}')

        # LOX.HP = -303200.254*1000, chamber_pressure
        # JetA.HP = -12978.768*1000, chamber_pressure

        print(f'gas init temp: {gas_init_temp}')

        self.gas.TP = gas_init_temp, chamber_pressure
        self.gas.equilibrate('HP')

        LOX.TP = LOX_inj_temp, LOX_inj_p
        LOX.equilibrate('HP')

        JetA.TP = JetA_inj_temp, JetA_inj_p
        JetA.equilibrate('HP')

        print(f'LOX: {LOX.TP}, JetA: {JetA.TP}, gas: {self.gas.TP}')

        # self.gas.HP = self.gas.enthalpy_mass, chamber_pressure
        # LOX.HP = LOX.enthalpy_mass, chamber_pressure
        # JetA.HP = JetA.enthalpy_mass, chamber_pressure

        self.gas.TP = self.gas.T, chamber_pressure
        LOX.TP = LOX.T, chamber_pressure
        JetA.TP = JetA.T, chamber_pressure

        print(f'LOX: {LOX.TP}, JetA: {JetA.TP}, gas: {self.gas.TP}')

        self.mixture.equilibrate('HP')

        print(f'LOX: {LOX.TP}, JetA: {JetA.TP}, gas: {self.gas.TP}')

        # chamber_pressure is 1.7976e+06 Pa

        # cp should be 3.3668666666666667
        # gamma should be 1.1849
        # temp should be 3041

        # LOX.SP = LOX.s, chamber_pressure
        # JetA.SP = JetA.s, chamber_pressure
        # self.gas.SP = self.gas.s, chamber_pressure
        # self.mixture.equilibrate('SP')
        gamma, cp = CanteraProperties._compute_gamma_cp(self.gas)

        print(f'cantera cp: {cp}, gamma: {gamma}')
        print(f'more direct: {self.gas.cp_mass}, gamma: {self.gas.cp_mass/self.gas.cv_mass}')
        self.mixture()

        error

        return self.gas.enthalpy_mass, gamma, entropy
    
    # https://kyleniemeyer.github.io/rocket-propulsion/thermochemistry/cea_cantera.html#rocket-calculations
    def get_throat_props(self,chamber_pressure,h_chamber):
        
        # Define Mach number vs. gamma equation to be solve numerically
        def iterate(pressure):
            self.gas.SP = self.chamber_entropy, pressure
            self.mixture.equilibrate('SP')
            gamma_throat, cp_throat = CanteraProperties._compute_gamma_cp(self.gas)
            speed_throat = np.sqrt(2*(h_chamber-cp_throat))
            Mach_throat = speed_throat / self.gas.sound_speed
            converged = np.abs(1-1/(Mach_throat**2)) < .4e-4
            next_pressure = self.gas.P*(1+gamma_throat*Mach_throat**2)/(1+gamma_throat)
            return converged, next_pressure, gamma_throat

        converged = False
        i = 0
        pressure_throat = chamber_pressure / np.power((self.chamber_gamma + 1) / 2., self.chamber_gamma / (self.chamber_gamma - 1))
        while not converged and i < self.max_iters_throat:
            converged, pressure_throat, gamma_throat = iterate(pressure_throat)
            i += 1

        return pressure_throat, self.gas.mean_molecular_weight / 1000, gamma_throat
    
    # Mostly copied & pasted from:
    # https://kyleniemeyer.github.io/rocket-propulsion/thermochemistry/cea_cantera.html#rocket-calculations
    @staticmethod
    def _compute_gamma_cp(gas):

        def get_thermo_derivatives(gas):
            '''Gets thermo derivatives based on shifting equilibrium.
            '''
            # unknowns for system with no condensed species:
            # dpi_i_dlogT_P (# elements)
            # dlogn_dlogT_P
            # dpi_i_dlogP_T (# elements)
            # dlogn_dlogP_T
            # total unknowns: 2*n_elements + 2

            num_var = 2 * gas.n_elements + 2

            coeff_matrix = np.zeros((num_var, num_var))
            right_hand_side = np.zeros(num_var)

            tot_moles = 1.0 / gas.mean_molecular_weight
            moles = gas.X * tot_moles

            # indices
            idx_dpi_dlogT_P = 0
            idx_dlogn_dlogT_P = idx_dpi_dlogT_P + gas.n_elements
            idx_dpi_dlogP_T = idx_dlogn_dlogT_P + 1
            idx_dlogn_dlogP_T = idx_dpi_dlogP_T + gas.n_elements

            # construct matrix of elemental stoichiometric coefficients
            stoich_coeffs = np.zeros((gas.n_elements, gas.n_species))
            for i, elem in enumerate(gas.element_names):
                for j, sp in enumerate(gas.species_names):
                    stoich_coeffs[i,j] = gas.n_atoms(sp, elem)

            # equations for derivatives with respect to temperature
            # first n_elements equations
            for k in range(gas.n_elements):
                for i in range(gas.n_elements):
                    coeff_matrix[k,i] = np.sum(stoich_coeffs[k,:] * stoich_coeffs[i,:] * moles)
                coeff_matrix[k, gas.n_elements] = np.sum(stoich_coeffs[k,:] * moles)
                right_hand_side[k] = -np.sum(stoich_coeffs[k,:] * moles * gas.standard_enthalpies_RT)

            # skip equation relevant to condensed species

            for i in range(gas.n_elements):
                coeff_matrix[gas.n_elements, i] = np.sum(stoich_coeffs[i, :] * moles)
            right_hand_side[gas.n_elements] = -np.sum(moles * gas.standard_enthalpies_RT)

            # equations for derivatives with respect to pressure

            for k in range(gas.n_elements):
                for i in range(gas.n_elements):
                    coeff_matrix[gas.n_elements+1+k,gas.n_elements+1+i] = np.sum(stoich_coeffs[k,:] * stoich_coeffs[i,:] * moles)
                coeff_matrix[gas.n_elements+1+k, 2*gas.n_elements+1] = np.sum(stoich_coeffs[k,:] * moles)
                right_hand_side[gas.n_elements+1+k] = np.sum(stoich_coeffs[k,:] * moles)

            for i in range(gas.n_elements):
                coeff_matrix[2*gas.n_elements+1, gas.n_elements+1+i] = np.sum(stoich_coeffs[i, :] * moles)
            right_hand_side[2*gas.n_elements+1] = np.sum(moles)
            
            derivs = np.linalg.solve(coeff_matrix, right_hand_side)

            dpi_dlogT_P = derivs[idx_dpi_dlogT_P : idx_dpi_dlogT_P + gas.n_elements]
            dlogn_dlogT_P = derivs[idx_dlogn_dlogT_P]
            dlogn_dlogP_T = derivs[idx_dlogn_dlogP_T]
            
            return dpi_dlogT_P, dlogn_dlogT_P, dlogn_dlogP_T

        dpi_dlogT_P, dlogn_dlogT_P, dlogn_dlogP_T = get_thermo_derivatives(gas)

        tot_moles = 1.0 / gas.mean_molecular_weight
        moles = gas.X * tot_moles
        
        # construct matrix of elemental stoichiometric coefficients
        stoich_coeffs = np.zeros((gas.n_elements, gas.n_species))
        for i, elem in enumerate(gas.element_names):
            for j, sp in enumerate(gas.species_names):
                stoich_coeffs[i,j] = gas.n_atoms(sp, elem)
        
        spec_heat_p = ct.gas_constant * (
            np.sum([dpi_dlogT_P[i] * 
                    np.sum(stoich_coeffs[i,:] * moles * gas.standard_enthalpies_RT) 
                    for i in range(gas.n_elements)
                    ]) +
            np.sum(moles * gas.standard_enthalpies_RT) * dlogn_dlogT_P +
            np.sum(moles * gas.standard_cp_R) +
            np.sum(moles * gas.standard_enthalpies_RT**2)
            )
        
        dlogV_dlogT_P = 1 + dlogn_dlogT_P
        dlogV_dlogP_T = -1 + dlogn_dlogP_T
        
        spec_heat_v = (
            spec_heat_p + gas.P * gas.v / gas.T * dlogV_dlogT_P**2 / dlogV_dlogP_T
            )

        gamma = spec_heat_p / spec_heat_v
        gamma_s = -gamma/dlogV_dlogP_T
        
        return gamma_s, spec_heat_p   

cant = CanteraProperties(P_inj=250*0.0689476*10**5,P_amb=1e5,con_rat=5.5,LOX_inj_temp=90.18,JetA_inj_temp=298.15,gas_init_temp=3000,MR_weight=1.8)

# full_species = {S.name: S for S in ct.Species.list_from_file('nasa_gas.yaml')}
# for species in full_species:
#     print(species)