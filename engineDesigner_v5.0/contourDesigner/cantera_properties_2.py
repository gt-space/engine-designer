import cantera as ct
from cantera.composite import SolutionArray
import numpy as np
import scipy.optimize as sci
import matplotlib.pyplot as plt

class Other:

    def __init__(self,mdot,P_inj,P_amb,con_rat, LOX_inj_temp, JetA_inj_temp,gas_init_temp,MR_weight,total_engine_length,max_iters_throat=1_000,max_iters=100):

        self.mdot = mdot
        self.P_inj = P_inj
        self.P_amb = P_amb
        self.con_rat = con_rat
        self.MR = MR_weight * 166 / 32 # ratio of moles
        self.total_engine_length = total_engine_length
        self.max_iters_throat = max_iters_throat
        self.max_iters = max_iters

        self.contour_defined = False

        # The below fields are initialized when set_contour is called
        self.area_arr = []
        self.throat_area = []
        self.pressure_arr = []
        self.temp_arr = []
        self.rho_arr = []
        self.ivac_arr = []
        self.cf_arr = []
        self.isp_arr = []
        self.h_arr = []
        self.e_arr = []
        self.molar_weight_arr = []
        self.cp_arr = []
        self.gamma_arr = []
        self.sound_speed_arr = []
        self.viscosity_arr = []
        self.thermal_conductivity_arr = []
        self.prandtl_arr = []
        self.velocity_arr = []

        all_species = {S.name: S for S in ct.Species.list_from_file('CombProducts.yaml',section='species')}
        self.JetA_X = {
            "C10H22(N)":  0.85,   # n-decane
            "C3H8":     0.15    # propane
        }

        self.LOX = ct.Solution(thermo='ideal-gas',species=[all_species['O2']])
        self.JetA = ct.Solution(thermo='ideal-gas',species=[all_species[S] for S in self.JetA_X])

        self.chamber_gas = ct.Solution(thermo='ideal-gas',kinetics='gas',transport_model='mixture-averaged',species=list(all_species.values()))
        self.throat = ct.Solution(thermo='ideal-gas',kinetics='gas',transport_model='mixture-averaged',species=list(all_species.values()))
        self.exit = ct.Solution(thermo='ideal-gas',kinetics='gas',transport_model='mixture-averaged',species=list(all_species.values()))
        reactions = ct.Reaction.list_from_file('CombProducts.yaml',self.chamber_gas)
        for reaction in reactions:
            self.chamber_gas.add_reaction(reaction)
            self.throat.add_reaction(reaction)
            self.exit.add_reaction(reaction)
        chamber_X_dict = {S: 0 for S in all_species}
        chamber_X_dict['O2'] = self.MR/(1+self.MR)
        for spec_name in self.JetA_X:
            chamber_X_dict[spec_name] = self.JetA_X[spec_name]/(1+self.MR)
        self.chamber_gas.X = chamber_X_dict

        self.mixed_propellants = ct.Mixture([(self.LOX,self.MR/(1+self.MR)),(self.JetA,1/(1+self.MR))])
        
        # https://rocketcea.readthedocs.io/en/latest/finite_area_comb.html
        self.Pc = P_inj/(1+.54/(con_rat**2.2))

        self._set_chamber_props(LOX_inj_temp,JetA_inj_temp)
        self._set_throat_props()

        # set expansion ratio
        k = self.throat_gamma # assume "frozen" at throat for expansion ratio calculation (reaction stops at the throat)
        self.exit_Mach = np.sqrt((2/(k-1))*(self.Pc/P_amb)**((k-1)/k)-1)
        self.exp_rat = (1/self.exit_Mach)*((1+(k+1)/2*self.exit_Mach**2)/((k+1)/2))**((k+1)/(2*(k-1)))

    @property
    def throat_PcOvPe(self):
        return self.throat_pressure / self.P_amb

    @property
    def get_eps(self):
        return self.exp_rat

    @property
    def get_exit_Mach(self):
        return self.exit_Mach

    @property
    def get_throat_MWt(self):
        return self.throat_MWt / 1000

    @property
    def get_throat_gamma(self):
        return self.throat_gamma

    @property
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
        self.chamber_area = self.area_arr[0]
        self.chamber_speed = self.mdot / (self.chamber_area*self.chamber_gas.density)
        i = 0
        while np.abs(self.area_arr[i] - self.chamber_area)/self.chamber_area<1e-4:
            i += 1
        self.chamber_end_ind = i -1
        i = 0
        while i < len(self.area_arr - 1) and self.area_arr[i] >= self.area_arr[i+1]:
            i += 1
        self.throat_ind = i
        self.throat_area = self.area_arr[self.throat_ind]
        p_conv = self._iterate_conv()
        self._set_props(p_conv)
        self.contour_defined = True
        print('CONTOUR SET')

    @property
    def pressures(self):
        if not self.area_arr:
            raise Exception('Contour must be initialized before property pressure is determined')
        return self.pressure_arr

    @property
    def temperatures(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property temperature is determined')
        return self.temp_arr
       
    @property 
    def pinj_ov_p(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property PinjOvP is determined')
        return self.P_inj / self.pressure_arr

    @property
    def Ae_ov_At(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property AeOvAt is determined')
        return self.area_arr / self.throat_area
        
    @property
    def cf(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property CF is determined')
        return self.cf_arr          

    @property
    def ivac(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property IVAC is determined')
        return self.ivac_arr

    @property
    def Isp(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property ISp is determined')
        return self.isp_arr
       
    @property 
    def rho(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property density is determined')
        return self.rho_arr        

    @property
    def h(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property enthalpy is determined')
        return self.h_arr

    @property
    def e(self): # return vector (internal energy)
        if not self.area_arr:
            raise Exception('Contour must be initialized before property internal energy is determined')
        return self.e_arr
        
    @property
    def MWt(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property molar weight is determined')
        return self.molar_weight_arr

    @property
    def cp(self): # return vector (specific heat capacity)
        if not self.area_arr:
            raise Exception('Contour must be initialized before property cp is determined')
        return self.cp_arr

    @property
    def gammas(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property gamma is determined')
        return self.gamma_arr
    
    @property
    def sound_speeds(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property sound speed is determined')
        return self.sound_speed_arr

    @property
    def viscosities(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property viscosity is determined')
        return self.viscosity_arr

    @property
    def conductivities(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property thermal conductivity is determined')
        return self.thermal_conductivity_arr

    @property
    def prandtl(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property prandtl number is determined')
        return self.prandtl_arr

    @property
    def velocity_contour(self): # return vector
        if not self.area_arr:
            raise Exception('Contour must be initialized before property velocity is determined')
        return self.velocity_arr
    
    def _set_chamber_props(self,LOX_inj_temp, JetA_inj_temp):

        '''
        CEA Output:
        REACTANT                    WT FRACTION                  ENERGY      TEMP
                                    (SEE note)                  KJ/KG-MOL      K
        FUEL        Jet-A(L)                     1.0000000   -303200.254    298.150
        OXIDANT     O2(L)                        1.0000000    -12978.768     90.180
        '''
        # Set initial propellant properties
        self.LOX.TP = LOX_inj_temp, self.Pc
        self.JetA.TP = JetA_inj_temp, self.Pc

        self.LOX()
        self.JetA()

        # Find mixed properties of combustion products before combustion
        self.mixed_propellants.equilibrate('HP')
        self.LOX()
        self.JetA()
        print(f'pc: {self.Pc:.2E}')

        # Set 1D reacting flame with inlet as end combustion
        self.chamber_gas.TP = self.JetA.TP[0], self.Pc
        print(self.chamber_gas.TP)
    #   #  self.chamber_gas.TP = 3e3, 20e5
    #     self.flame = ct.FreeFlame(gas=self.chamber_gas, width=self.total_engine_length)
      
    #     init_guess = SolutionArray(self.chamber_gas)
    # #    init_guess.__setattr__('t',[self.chamber_gas.T, self.chamber_gas.T*1.2, self.chamber_gas.T*2])
    # #    init_guess.__setattr__('T', self.chamber_gas.T*1.2) # TODO : Make this a ternary w / optional user input
    #     self.flame.set_initial_guess(locs=[0, .25, .5, .75, 1], data=init_guess)
    
    #     self.flame.solve()
        self.chamber_gas.equilibrate('HP')
        gamma, _ = Other._compute_gamma_cp(self.chamber_gas)
        self.chamber_rho = self.chamber_gas.density
        self.chamber_temp = self.chamber_gas.T

        self.h_chamber = self.chamber_gas.enthalpy_mass
        self.chamber_gamma = gamma
        self.chamber_entropy = self.chamber_gas.entropy_mass
    
    # https://kyleniemeyer.github.io/rocket-propulsion/thermochemistry/cea_cantera.html#rocket-calculations
    def _set_throat_props(self):
        self.throat.TPX = self.chamber_gas.T, self.chamber_gas.P, self.chamber_gas.X
        self.throat.SPX = self.chamber_gas.entropy_mass, self.chamber_gas.P, self.chamber_gas.X
        
        # Define Mach number vs. gamma equation to be solve numerically
        def iterate(pressure):
            self.throat.SP = self.chamber_entropy, pressure
            self.throat.equilibrate('SP')
            gamma_throat, cp_throat = Other._compute_gamma_cp(self.throat)
            speed_throat = np.sqrt(2*(self.h_chamber-self.throat.enthalpy_mass))
            Mach_throat = speed_throat / self.throat.sound_speed
            converged = np.abs(1-1/(Mach_throat**2)) < .4e-4
            next_pressure = self.throat.P*(1+gamma_throat*Mach_throat**2)/(1+gamma_throat)
            return converged, next_pressure, gamma_throat, cp_throat

        converged = False
        i = 0
        pressure_throat = self.Pc / np.power((self.chamber_gamma + 1) / 2., self.chamber_gamma / (self.chamber_gamma - 1))
        print(f'pressure throat: {pressure_throat:.2E} vs p amb: {self.P_amb:.2E} ')
        while not converged and i < self.max_iters_throat:
            converged, pressure_throat, gamma_throat, cp_throat = iterate(pressure_throat)
            i += 1
        
        if not converged:
            print(f'Warning: solver did not converge and throat Mach number is {np.sqrt(2*(self.h_chamber-self.throat.enthalpy_mass))/self.throat.sound_speed}')

        self.rho_throat = self.throat.density
        self.throat_pressure = pressure_throat
        self.throat_MWt = self.throat.mean_molecular_weight
        self.throat_gamma = gamma_throat
        self.throat_temp = self.throat.T
        self.h_throat = self.throat.enthalpy_mass
        self.e_throat = self.throat.int_energy_mass
        self.cp_throat = cp_throat 
        self.cv_throat = self.cp_throat/self.throat_gamma   

    # Find the pressure at point of the converging section using Cantera equilibrium solver  
    def _iterate_conv(self):
        self.conv_areas = self.area_arr[self.chamber_end_ind:self.throat_ind+1]
        self.div_areas = self.area_arr[self.throat_ind+1:]
        self.exit.TPX = self.chamber_gas.T, self.chamber_gas.P, self.chamber_gas.X
        self.exit.SPX = self.chamber_gas.entropy_mass, self.chamber_gas.P, self.chamber_gas.X
        self.exit.SP = self.chamber_entropy, self.P_amb
        self.exit.equilibrate('SP')
        p_conv = np.linspace(self.Pc,self.throat_pressure,self.throat_ind-self.chamber_end_ind+1)
        func = lambda p : np.sum(np.square(self._p_to_A(p) - self.conv_areas))
        bounds = sci.Bounds(self.P_amb/self.Pc,1)
        sol = sci.minimize(func, p_conv/self.Pc,bounds=bounds)
        print(sol.message)
        return np.array(sol.x) * self.Pc

    # Given a normalized pressure array, find the area of the converging section at each point
    def _p_to_A(self, p_arr_norm):
        p_arr = p_arr_norm*self.Pc
        u = self.chamber_speed
        gas = self.chamber_gas
        rho_arr = []
        A_arr = [self.chamber_area]
        for i in range(len(p_arr_norm)):
            gas.SP = self.chamber_entropy, p_arr[i]
            gas.equilibrate('SP')
            rho_arr.append(gas.density)
        for i in range(1, len(p_arr_norm)):
            rho = rho_arr[i]
            if rho_arr[i] != rho_arr[i-1]:
                u_du = -(p_arr[i]-p_arr[i-1]) # momentum transport
                u_ov_du = -rho_arr[i]/(rho_arr[i]-rho_arr[i-1]) # continuity
                du = np.sqrt(u_du/u_ov_du)
                u += du
                A = self.mdot/(rho*u)
                A_arr.append(A)
              #  print(f'r: {np.sqrt((2/np.pi)*A)/.0254:.2E}, dA: {A-self.nozzle_areas[i-1]/.0254**2:.2E}')
            else:
                A_arr.append(A_arr[i-1])
        return A_arr

    def _set_props(self, conv_pressure):
        '''
        Functions for determining ivac (specific impulse in vaccum), cf (thrust coefficient), isp (specific impulse),
        gamma, viscosity, thermal conductivity, and prandtal number
        '''
        gas = self.chamber_gas

        def compute_ivac(Isp,temp,R,Mwt):
            return Isp + temp * R / (Isp * Mwt)
        
        def compute_cstar(gamma, temperature, molecular_weight):
            return (
                np.sqrt(ct.gas_constant * temperature / (molecular_weight * gamma)) *
                np.power(2 / (gamma + 1), -(gamma + 1) / (2*(gamma - 1)))
                )

        def compute_isp(enthalpy):
            return np.sqrt(2 * (self.h_chamber - enthalpy))

        def compute_viscosity():
            pass

        def compute_thermal_conductivity():
            pass

        def compute_prandtl():
            pass

        # Converging section
        for i in range(len(conv_pressure)):
            p = conv_pressure[i]
            gas.SP = self.chamber_entropy, p
            gamma, cp = Other._compute_gamma_cp(gas)
            a = np.sqrt(gamma*ct.gas_constant*gas.T/gas.mean_molecular_weight)
            Isp = compute_isp(gas.enthalpy_mass)
            self.temp_arr.append(gas.T)
            self.rho_arr.append(gas.density_mass)
            self.ivac_arr.append(compute_ivac(Isp,gas.T,ct.gas_constant/gas.mean_molecular_weight,gas.mean_molecular_weight))
            self.cf_arr.append(Isp / compute_cstar(gamma,gas.T,gas.mean_molecular_weight/1000))
            self.isp_arr.append(Isp)
            self.h_arr.append(gas.enthalpy_mass)
            self.e_arr.append(gas.int_energy_mass)
            self.molar_weight_arr.append(gas.mean_molecular_weight/1000) # kg/kmol to kg/mol
            self.cp_arr.append(cp)
            self.gamma_arr.append(gamma)
            self.sound_speed_arr.append(a)
            self.velocity_arr.append(self.mdot/(gas.density_mass*self.conv_areas[i]))
            self.viscosity_arr.append(0) # TBD
            self.thermal_conductivity_arr.append(0) # TBD
            self.prandtl_arr.append(0) # TBD
        self.pressure_arr=conv_pressure.tolist()

        def Mach_func(Mach):
            exp=(self.throat_gamma+1)/(2*(self.throat_gamma-1))
            t1=(self.throat_gamma+1)/2
            t2=1+Mach**2*(self.throat_gamma-1)/2
            return t1**(-exp)*(t2**exp)/Mach

        # Diverging section: assume frozen flow, CPG
        for A in self.div_areas:
            Mach = sci.root(lambda Mach: Mach_func(Mach)-A/self.throat_area, self.exit_Mach).x[0]
            print(f'Mach: {Mach}, area ratio: {A/self.throat_area}')
            T=self.throat_temp*(1+Mach**2*(self.throat_gamma-1)/2)**(-1)
            p=self.throat_pressure*(T/self.throat_temp)**(self.throat_gamma/(self.throat_gamma-1))
            rho=self.rho_throat*(p/self.throat_pressure)**(1/self.throat_gamma)
            h=self.h_throat+(T-self.throat_temp)*self.cp_throat
            Isp = compute_isp(h)
            self.pressure_arr.append(p)
            self.temp_arr.append(T)
            self.rho_arr.append(rho)
            self.ivac_arr.append(compute_ivac(Isp,T,ct.gas_constant/self.throat_MWt,self.throat_MWt))
            self.cf_arr.append(Isp / compute_cstar(self.throat_gamma,T,self.throat_MWt/1000))
            self.isp_arr.append(Isp)
            self.h_arr.append(h)
            self.e_arr.append(self.e_throat+self.cv_throat*(T-self.throat_temp))
            self.molar_weight_arr.append(self.throat_MWt/1000) # kg/kmol to kg/mol
            self.cp_arr.append(cp)
            self.gamma_arr.append(self.throat_gamma)
            self.sound_speed_arr.append(np.sqrt(self.throat_gamma*ct.gas_constant*T/self.throat_MWt))
            self.velocity_arr.append(Mach*self.throat.sound_speed)
            self.viscosity_arr.append(0) # TBD
            self.thermal_conductivity_arr.append(0) # TBD
            self.prandtl_arr.append(0) # TBD

        # Add chamber properties and turn everything into a numpy array
        self.pressure_arr = np.concat([np.ones(self.chamber_end_ind)*self.Pc,np.array(self.pressure_arr)])
        self.temp_arr = np.concat([np.ones(self.chamber_end_ind)*self.chamber_temp,np.array(self.temp_arr)])
        self.rho_arr = np.concat([np.ones(self.chamber_end_ind)*self.chamber_rho,np.array(self.rho_arr)])
        self.ivac_arr = np.concat([np.ones(self.chamber_end_ind)*self.ivac_arr[0],np.array(self.ivac_arr)])
        self.cf_arr = np.concat([np.ones(self.chamber_end_ind)*self.cf_arr[0],np.array(self.cf_arr)])
        self.isp_arr = np.concat([np.ones(self.chamber_end_ind)*self.isp_arr[0],np.array(self.isp_arr)])
        self.h_arr = np.concat([np.ones(self.chamber_end_ind)*self.h_chamber,np.array(self.h_arr)])
        self.e_arr = np.concat([np.ones(self.chamber_end_ind)*self.e_arr[0],np.array(self.e_arr)])
        self.molar_weight_arr = np.concat([np.ones(self.chamber_end_ind)*self.molar_weight_arr[0],np.array(self.molar_weight_arr)])
        self.cp_arr = np.concat([np.ones(self.chamber_end_ind)*self.cp_arr[0],np.array(self.cp_arr)])
        self.gamma_arr = np.concat([np.ones(self.chamber_end_ind)*self.chamber_gamma,np.array(self.gamma_arr)])
        self.sound_speed_arr = np.concat([np.ones(self.chamber_end_ind)*self.sound_speed_arr[0],np.array(self.sound_speed_arr)])
        self.viscosity_arr = np.concat([np.ones(self.chamber_end_ind)*self.viscosity_arr[0],np.array(self.viscosity_arr)])
        self.thermal_conductivity_arr = np.concat([np.ones(self.chamber_end_ind)*self.thermal_conductivity_arr[0],np.array(self.thermal_conductivity_arr)])
        self.prandtl_arr = np.concat([np.ones(self.chamber_end_ind)*self.prandtl_arr[0],np.array(self.prandtl_arr)])
        self.velocity_arr = np.concat([np.ones(self.chamber_end_ind)*self.chamber_speed,np.array(self.velocity_arr)])
        self.Mach_arr = self.velocity_arr / self.sound_speed_arr
    
    
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
    
cant = Other(mdot=20,P_inj=250*0.0689476*10**5,P_amb=1e5,con_rat=5.5,\
                         LOX_inj_temp=90.18,JetA_inj_temp=298.15,\
                            gas_init_temp=3000,MR_weight=1.8,total_engine_length=1)

print(f'exit mach: {cant.exit_Mach}')
chm = 10
thr = 4
print(f'val: {thr*cant.exp_rat}')
A=.0254*np.concat([chm*np.ones(50),np.linspace(chm,thr,25),np.linspace(thr,thr*5,25)])
cant.set_contour(A)

x=np.linspace(0,1,len(A))

fig,ax=plt.subplots(4,1)
ax[0].plot(x/.0254, A/.0254**2)
ax[0].set_xlabel('Axial Position (inches)')
ax[1].set_xlabel('Axial Position (inches)')
ax[2].set_xlabel('Axial Position (inches)')
ax[3].set_xlabel('Axial Position (inches)')
ax[1].plot(x/.0254, cant.Mach_arr)
ax[2].plot(x/.0254,cant.pressure_arr)
ax[3].plot(x/.0254,cant.temp_arr)
ax[0].set_title('Area (in^2)')
ax[1].set_title('Mach Number')
ax[2].set_title('Pressure (Pa)')
ax[3].set_title('Temperature (K)')
plt.show()