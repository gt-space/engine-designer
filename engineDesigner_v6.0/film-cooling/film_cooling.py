# main.py
import json
import thermo_properties as thermo
import hydrodynamics as hydro
import heat_mass_transfer as hmt
import post_processing as post

# 1. Initialization
config = json.load("config.json")
x = 0.0
dx = config.dx_m
Gamma_L = config.Gamma_L_initial
T_L = config.T_L_initial

# Data storage for plotting
data_log = []

# 2. Main Space-Marching Loop
while Gamma_L > 0 and x < config.chamber_length_m:
    
    # Get Properties
    liquid_props = thermo.get_liquid_props(T_L, config.P_cc_bar)
    gas_props = thermo.get_gas_props(config.T_g_K, config.P_cc_bar)
    P_D = 0.5 * gas_props.rho_g * (config.U_g_m_s ** 2) #Dynamic Pressure
    C_f = hydro.calculate_friction_factor(gas_props)
    
    # Hydrodynamics
    delta_L, U_ls = hydro.solve_couette(Gamma_L, liquid_props, P_D, C_f)
    Fr, freq, wave_len = hydro.check_stability_and_roughness(config.U_g_m_s, U_ls, delta_L, gas_props, liquid_props, config.G_g)
    
    # Heat and Mass Transfer
    m_dot_Ent = hmt.calc_entrainment(P_D, config.U_g_m_s, liquid_props, delta_L, C_f)
    q_Tot, m_dot_Evap = hmt.calc_heat_and_evaporation(T_L, Fr, gas_props, liquid_props, config, C_f)
    
    # 3. ODE Integration (Euler Method)
    if T_L < liquid_props.T_sat:
        # Energy goes into raising temperature
        dT_dx = q_Tot / (Gamma_L * liquid_props.Cp_L)
        T_L = T_L + (dT_dx * dx)
        # need to add some correction line if T_L overshoots saturation temp.
    elif T_L >= liquid_props.T_sat: 
            T_L = liquid_props.T_sat
    
    # Subtract mass (both entrainment and evaporation)
    dGamma_dx = - (m_dot_Ent + m_dot_Evap)
    Gamma_L = Gamma_L + (dGamma_dx * dx)
    
    # Save step data
    data_log.append({x, T_L, Gamma_L, delta_L, m_dot_Ent, m_dot_Evap})
    
    # Step forward
    x = x + dx

# 4. End Simulation
print("Coolant Depleted!")
print("Film-Cooled Length (FCL) = ", x, " meters")

post.plot_results(data_log)