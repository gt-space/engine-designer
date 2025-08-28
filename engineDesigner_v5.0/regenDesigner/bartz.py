import numpy as np
import sys

sys.path.insert(0, '/Users/atand/OneDrive/Documents/Code stuff/engine-designer/engineDesigner_v5.0/regenDesigner')

from contourDesigner.CEA_properties import ceaToSI
from utils.fuel_props import JetA

def bartz(engine, T_wg, i):
    # DUE TO BARTZ CORRELATION BEING DEVELOPED IN ENGLISH UNITS, VALUES
    # CONVERTED TO ENGLISH UNITS IN EQUATIONS BELOW THEN RECONVERTED

    recovery_factor = 1
    cstar_eff = engine.cstar_eff

    g = 9.8104 # Gravity, m/s^2
    C_star = engine.C_star
    R_tCurve = engine.R_tCurve
    R_t = engine.R_t
    P_c = engine.P_inj_psi * 6894.76 # Injector face pressure in psi
    T_c0 = engine.engineProps[0,9] # Convert K to R
    gam0 = engine.engineProps[0,15] # Chamber stagnation gamma
    # Average of frozen and not frozen
    cp_ns = (engine.engineProps[99,20] + engine.engineProps[99,14]) * 0.5 * 1000#Convert from Joules to BTU/lb*F
    praneff_ns = 0.5 * (engine.engineProps[99,22] + engine.engineProps[99,19]) # No conversion needed
    '''
    praneff_ns_f = engine.engineProps[99,22] # No conversion needed
    praneff_ns_eq = engine.engineProps[99,19]

    visc_ns = (engine.engineProps[0,17]/1000) * 0.1 # Dyn viscosity (mili-Poise > Poise > kg/ms)
    cond = 0.1 * (engine.engineProps[i,21] + engine.engineProps[i,18]) * 0.5 # Eq Conductivity (W/mK)
    cp_ns_calc = (engine.engineProps[99,20] + engine.engineProps[99,14]) * 1000 #(kJ > J)
    # praneff_ns = visc_ns*cp_ns_calc/cond
    # praneff_ns = 0.508
    # praneff_ns_est = (4*gam0)/(9*gam0-5) #Originally used effective prandtl number from CEA, switched to gamma correlation
    '''
    visc_ns = (engine.engineProps[0,17]/1000) * 0.1# Convert to Poise

    # Run Bartz Correlation
    gam = engine.engineProps[i, 15] #Index gamma; No conversion needed
    mach = engine.engineProps[i, 4]
    contourR = engine.engineProps[i, 0] #Convert to Inches
    T_aw = T_c0 * (1 + (praneff_ns**(1/3)) * ((gam - 1)/2) * (mach**2))/(1 + ((gam - 1)/2) * (mach**2))
    T_aw = T_aw * cstar_eff**2
    #Sigma = sigA * sigB - Split into Terms for readability
    sigA = (0.5 * (T_wg/T_c0) * (1 + ((gam - 1)/2) * mach**2) + 0.5)**(-0.68)
    sigB = (1 + ((gam - 1)/2) * (mach**2))**(-0.12)
    sigma = sigA*sigB
    # a_T = (aA*aB*aC*aD)*aE*sigma - Split into Terms
    aA = 0.026/((2 * R_t)**0.2)
    aB = ((visc_ns**0.2) * cp_ns)/(praneff_ns**0.6)
    aC = (P_c/C_star)**0.8
    aD = ((2 * R_t)/R_tCurve)**0.1
    aE = ((R_t**2)/(contourR**2))**0.9
    h_g = (aA * aB * aC * aD) * aE * sigma * 0.75 #Convection coefficient with soot knockdown
    q_conv = h_g * (T_aw - T_wg) #Heat flux (heat per area)

    # Adiabtic wall temp T_aw is temp if the wall was adiabatic (no heat transfer
    # this means that this would be equal to the flow temperature at that region
    # thus, T_aw = T_hg

    # print(h_g)
    return (h_g, q_conv, T_aw)

def hg_gas_film(engine, T_wg, T_hg, deltaQ, last_wall_temp, last_u_cool, M_wt, i, dz):
    # T_wg assumed to be coolant temperature at wall, T_hg is mainstream gas temperature
    # deltaQ is negative value for specific heat transfer rate out of coolant from last axial station
    # L is circumference (m)
    # This function is only applied when T_wg is not T_hg

    # print(f"T_wg: {T_wg}, T_hg: {T_hg}, last_u_cool: {last_u_cool}, M_wt: {M_wt}, dz: {dz}")

    X = i * dz # length downstream of coolant injection (m)

    # get gas properties 
    gamma_comb_products = engine.engineProps[i, 15]
    R = 8.314/M_wt # specific gas constant of combustion products (kJ/kgK)
    u_hg = engine.engineProps[i, 4]*np.sqrt(gamma_comb_products*T_hg*R)
    pressure = engine.engineProps[i, 8]*10**5 # bars to Pascals

    # update coolant properties
    cp_cool = JetA.get_cp_vapor(T_wg)
    alpha_cool = JetA.get_conductivity(T_wg)/(JetA.get_rho_v(T_wg, pressure)*cp_cool) # thermal diffusivity
    # deltaT_cool = T_wg-last_wall_temp
    # if i > 0:
    #     deltaT_cool = engine.engineProps[i,9]-engine.engineProps[i-1,9]
    #     print(f"deltaT_cool: {deltaT_cool}")
    # else:
    #     deltaT_cool = 0
    # deltaK = deltaQ + cp_cool * deltaT_cool # deltaQ is net heat flux into coolant & the second term is heat flux to mainstream gases; this solves for change in specific kinetic energy of film
    deltaK = 0 # should fix the above later
    u_cool = np.sqrt(2*(.5*last_u_cool**2+deltaK))

    h_g = engine.film_cooling[-1][0].get_h_g(u_hg, u_cool, cp_cool, alpha_cool, i, X)
    # print(f'hg: {h_g}, T_hg: {T_hg}, T_wg: {T_wg}, z: {i*dz}')

    return h_g, u_cool, T_hg > T_wg+5

def hg_boiling_liquid_film(engine, T_wl, T_hg, heat_flux_wall, M_wt, i, dz):
    # T_wl assumed to be liquid coolant temperature at wall, T_hg is mainstream gas temperature
    # heat_flux_wall is positive value for the current iteration's heat flux through the wall

    liquid_film = engine.film_cooling[-1][0] # instance of liquid_film_cooling
    radius = engine.engineProperties[i][0]*.3048 # ft to m
    mdot_cool = engine.film_cooling[1]

    # get hot gas properties
    gamma_comb_products = engine.engineProps[i, 15]
    R = 8.314/ceaToSI(M_wt, ) # specific gas constant of combustion products (kJ/kgK)
    u_hg = engine.engineProps[i, 4]*np.sqrt(gamma_comb_products*T_hg*R) # m/s
    pressure = engine.engineProps[i, 8]*10**5 # bars to Pascals
    mew_gas = ceaToSI(engine.engineProps[i,17])
    pr_gas = engine.engineProps[i, 19]
    cp_gas = ceaToSI(engine.engineProps[i, 14])

    # update liquid coolant properties
    saturation_temp = JetA.get_saturation_temp(pressure) # K
    h_fg = JetA.get_h_fg(saturation_temp) # J/kgK
    temp_avg_cool = (T_wl+saturation_temp)/2 # average temperature between boiling surface & wall
    rho_v_cool = JetA.get_rho_v(T_hg) # coolant vapor density, kg/m^3
    c_l_cool = JetA.get_c_liquid(temp_avg_cool)
    hstar_fg = h_fg + (saturation_temp-temp_avg_cool)*c_l_cool
    
    enthalpy_flux_other = -heat_flux_wall - T_wl*c_l_cool*mdot_cool/(2*np.pi*radius)*dz #see liquid_film_cooling.py for details

    h_g = liquid_film.get_h_g(1, hstar_fg, u_hg, rho_v_cool, mew_gas, M_wt, pr_gas, cp_gas, T_hg, saturation_temp, enthalpy_flux_other, liquid_film.radii(i)*2)
    
    return h_g


def bartz_sigma(engine, i, recovery_factor=1):
    #Computes bartz correlation without computing the boundary layer correction factor (sigma)
    #Used for solving for the inner wall temperature
    
    recovery_factor = 1
    cstar_eff = engine.cstar_eff

    C_star = engine.C_star
    R_tCurve = engine.R_tCurve
    R_t = engine.R_t
    P_c = engine.P_inj_psi * 6894.76# Injector face pressure in psi
    T_c0 = engine.engineProps[0,9]
    # gam0 = engine.engineProps[0,15] # Chamber stagnation gamma
    # Average of frozen and not frozen
    cp_ns = (engine.engineProps[99,20] + engine.engineProps[99,14]) * 0.5 * 1000
    praneff_ns = 0.5 * (engine.engineProps[99,22] + engine.engineProps[99,19]) # No conversion needed
    # praneff_ns_f = engine.engineProps[99,22] # No conversion needed
    # praneff_ns_eq = engine.engineProps[99,19]

    visc_ns = (engine.engineProps[0,17]/1000) * 0.1 # Dyn viscosity (mili-Poise > Poise > kg/ms)
    # cond = 0.1 * (engine.engineProps[i,21] + engine.engineProps[i,18]) * 0.5 # Eq Conductivity (W/mK)
    # cp_ns_calc = (engine.engineProps[99,20] + engine.engineProps[99,14]) * 1000 #(kJ > J)
    # praneff_ns = visc_ns*cp_ns_calc/cond
    # praneff_ns = 0.508
    # praneff_ns_est = (4*gam0)/(9*gam0-5) #Originally used effective prandtl number from CEA, switched to gamma correlation




    # Run Bartz Correlation
    gam = engine.engineProps[i, 15] #Index gamma; No conversion needed
    mach = engine.engineProps[i, 4]
    contourR = engine.engineProps[i, 0]
    T_aw = T_c0 * (1 + (praneff_ns**(1/3)) * ((gam - 1)/2) * (mach**2))/(1 + ((gam - 1)/2) * (mach**2))
    T_aw = T_aw * cstar_eff**2
    #Sigma = sigA * sigB - Split into Terms for readability

    aA = 0.026/((2 * R_t)**0.2)
    aB = ((visc_ns**0.2) * cp_ns)/(praneff_ns**0.6)
    aC = (P_c/C_star)**0.8
    aD = ((2 * R_t)/R_tCurve)**0.1
    aE = ((R_t**2)/(contourR**2))**0.9
    h_g = (aA * aB * aC * aD) * aE * 0.5 #Convection coefficient

    # print(h_g)
    return (h_g, T_aw)
