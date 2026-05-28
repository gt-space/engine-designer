# heat_mass_transfer.py

import math


STEFAN_BOLTZMANN = 5.670374419e-8  # W/(m^2·K^4)

def calc_entrainment(P_D, U_g, gas_props, liquid_props, delta_L, C_f):
    # Calculate Entrainment Parameter (X_e)
    T_ratio = gas_props.T_g / liquid_props.T_sat
    X_e = ((P_D ** 0.5) / liquid_props.sigma_L) * (T_ratio ** 0.25)
    
    if X_e >= 474.14: # 2025 SI Threshold
        # Calculate mass susceptible to tearing
        m_dot_cool = (P_D * C_f * (delta_L ** 2) * liquid_props.rho_L) / (2 * liquid_props.mu_L)
        # Calculate mass actually stripped
        m_dot_Ent = 1.41e-3 * (X_e - 474.14) * m_dot_cool
    else:
        m_dot_Ent = 0.0
        
    return m_dot_Ent

def calc_heat_and_evaporation(T_L, Fr, gas_props, liquid_props, config, C_f):
    # 1. Radiation
    q_Rad = STEFAN_BOLTZMANN * config.emissivity * ((gas_props.T_g ** 4) - (T_L ** 4)) #stefan-boltzman constant 
    
    # 2. Baseline Convection (with roughness penalty)
    St_0 = 0.5 * C_f * (gas_props.Pr_g ** -0.667) * (1 + Fr)
    K_Tu = 1 + (4 * config.turbulence_intensity_eT)
    h_0 = K_Tu * config.G_g * gas_props.Cp_g * St_0
    
    # 3. Phase check
    if T_L < liquid_props.T_sat:
        # PRE-SATURATION: Heating up, no boiling yet
        m_dot_Evap = 0.0
        q_Conv = h_0 * (gas_props.T_g - T_L)
        q_Tot = q_Rad + q_Conv
        
    else:
        # POST-SATURATION: Iterative solver for vapor blowing effect
        T_L = liquid_props.T_sat
        m_dot_Evap_guess = 0.1 # Initial guess
        error = 1.0
        
        while error > 0.001:
            # Calculate blowing ratio (H)
            H = gas_props.Cp_g * config.KM * (m_dot_Evap_guess / h_0)
            
            # Reduce heat transfer coefficient due to vapor shield
            h_blown = h_0 * math.ln((1 + H) / H)
            
            # Re-calculate heat fluxes
            q_Conv = h_blown * (gas_props.T_g - T_L)
            q_Tot_new = q_Rad + q_Conv
            
            # Re-calculate evaporation based on new heat flux
            m_dot_Evap_new = q_Tot_new / liquid_props.LHV
            
            # Check convergence
            error = abs(m_dot_Evap_new - m_dot_Evap_guess)
            m_dot_Evap_guess = m_dot_Evap_new
            
        m_dot_Evap = m_dot_Evap_new
        q_Tot = q_Tot_new

    return q_Tot, m_dot_Evap