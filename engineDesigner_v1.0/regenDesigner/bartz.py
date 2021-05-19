
def bartz(engine, T_wg, i):
    # DUE TO BARTZ CORRELATION BEING DEVELOPED IN ENGLISH UNITS, VALUES
    # CONVERTED TO ENGLISH UNITS IN EQUATIONS BELOW THEN RECONVERTED
    g = 32.174 # Gravity, ft/s^2
    T_wg = T_wg * 1.8 # Convert K to R (initial guess for finding coolant temps)
    C_star = engine.C_star * 3.28084 # Convert m/s to ft/s
    R_tCurve = engine.R_tCurve * 39.3701 # Convert m to in
    R_t = engine.R_t * 39.3701 # Convert m to in
    P_c = engine.P_inj_psi # Injector face pressure in psi
    T_c0 = engine.engineProps[0,9] * 1.8 # Convert K to R
    gam0 = engine.engineProps[0,15] # Chamber stagnation gamma
    # Average of frozen and not frozen
    cp_ns = (engine.engineProps[99,20] + engine.engineProps[99,14]) * 0.5 * 0.23884; #Convert from Joules to BTU/lb*F
    praneff_ns = 0.5 * (engine.engineProps[99,22] + engine.engineProps[99,19]) # No conversion needed
    praneff_ns_f = engine.engineProps[99,22] # No conversion needed
    praneff_ns_eq = engine.engineProps[99,19]

    visc_ns = (engine.engineProps[0,17]/1000) * 0.1 # Dyn viscosity (mili-Poise > Poise > kg/ms)
    cond = 0.1 * (engine.engineProps[i,21] + engine.engineProps[i,18]) * 0.5 # Eq Conductivity (W/mK)
    cp_ns_calc = (engine.engineProps[99,20] + engine.engineProps[99,14]) * 1000 #(kJ > J)
    # praneff_ns = visc_ns*cp_ns_calc/cond
    # praneff_ns = 0.508
    praneff_ns_est = (4*gam0)/(9*gam0-5) #Originally used effective prandtl number from CEA, switched to gamma correlation

    visc_ns = (engine.engineProps[0,17]/1000) * (0.0672/12) # Convert to Poise then to lb/in-s

    bartzData = [T_wg, T_c0, R_t, P_c, C_star, R_tCurve, praneff_ns, visc_ns, cp_ns, g]

    # unpack data
    T_wg = bartzData[0]
    T_c0 = bartzData[1]
    T_c0 = 3008 * 1.8
    R_t = bartzData[2]
    P_c = bartzData[3]
    C_star = bartzData[4]
    R_tCurve = bartzData[5]
    praneff_ns = bartzData[6]
    visc_ns = bartzData[7]
    cp_ns = bartzData[8]
    g = bartzData[9]

    # Run Bartz Correlation
    gam = engine.engineProps[i, 15] #Index gamma; No conversion needed
    mach = engine.engineProps[i, 5]
    contourR = engine.engineProps[i, 0] * 39.3701 #Convert to Inches
    T_aw = T_c0 * (1 + (praneff_ns**(1/3)) * ((gam - 1)/2) * (mach**2))/(1 + ((gam - 1)/2) * (mach**2))

    #Sigma = sigA * sigB - Split into Terms for readability
    sigA = (0.5 * (T_wg/T_c0) * (1 + ((gam - 1)/2) * mach**2) + 0.5)**(-0.68)
    sigB = (1 + ((gam - 1)/2) * (mach**2))**(-0.12)
    sigma = sigA*sigB
    # a_T = (aA*aB*aC*aD)*aE*sigma - Split into Terms
    aA = 0.026/((2 * R_t)**0.2)
    aB = ((visc_ns**0.2) * cp_ns)/(praneff_ns**0.6)
    aC = (P_c * g/C_star)**0.8;
    aD = ((2 * R_t)/R_tCurve)**0.1
    aE = ((R_t**2)/(contourR**2))**0.9
    h_g = (aA * aB * aC * aD) * aE * sigma #Convection coefficient
    q_conv = h_g * (T_aw - T_wg) #Heat flux (heat per area)

    # Re-Unit Convert
    h_g = h_g * 144 * 20428.175 #BTU/in^2-s-F -> BTU/ft^2-s-F -> W/m^2K
    q_conv = q_conv * 1634246.235 #BTU/in^2-s -> W/m^2
    T_aw = T_aw / 1.8 #R -> K

    # print(h_g)
    return (h_g, q_conv, T_aw)
    # Run Bartz Correlation
    gam = engine.engineProps[i, 15] #Index gamma; No conversion needed
    mach = engine.engineProps[i, 5]
    contourR = engine.engineProps[i, 0] * 39.3701 #Convert to Inches
    T_aw = T_c0 * (1 + (praneff_ns**(1/3)) * ((gam - 1)/2) * (mach**2))/(1 + ((gam - 1)/2) * (mach**2))

    #Sigma = sigA * sigB - Split into Terms for readability
    sigA = (0.5 * (T_wg/T_c0) * (1 + ((gam - 1)/2) * mach**2) + 0.5)**(-0.68)
    sigB = (1 + ((gam - 1)/2) * (mach**2))**(-0.12)
    sigma = sigA*sigB
    # a_T = (aA*aB*aC*aD)*aE*sigma - Split into Terms
    aA = 0.026/((2 * R_t)**0.2)
    aB = ((visc_ns**0.2) * cp_ns)/(praneff_ns**0.6)
    aC = (P_c * g/C_star)**0.8;
    aD = ((2 * R_t)/R_tCurve)**0.1
    aE = ((R_t**2)/(contourR**2))**0.9
    h_g = (aA * aB * aC * aD) * aE * sigma #Convection coefficient
    q_conv = h_g * (T_aw - T_wg) #Heat flux (heat per area)

    # Re-Unit Convert
    h_g = h_g * 144 * 20428.175 #BTU/in^2-s-F -> BTU/ft^2-s-F -> W/m^2K
    q_conv = q_conv * 1634246.235 #BTU/in^2-s -> W/m^2
    T_aw = T_aw / 1.8 #R -> K

    # print(h_g)
    return (h_g, q_conv, T_aw)
    T_aw = T_aw / 1.8 #R -> K

    # print(h_g)
    return (h_g, q_conv, T_aw)
    # Run Bartz Correlation
    gam = engine.engineProps[i, 15] #Index gamma; No conversion needed
    mach = engine.engineProps[i, 5]
    contourR = engine.engineProps[i, 0] * 39.3701 #Convert to Inches
    T_aw = T_c0 * (1 + (praneff_ns**(1/3)) * ((gam - 1)/2) * (mach**2))/(1 + ((gam - 1)/2) * (mach**2))

    #Sigma = sigA * sigB - Split into Terms for readability
    sigA = (0.5 * (T_wg/T_c0) * (1 + ((gam - 1)/2) * mach**2) + 0.5)**(-0.68)
    sigB = (1 + ((gam - 1)/2) * (mach**2))**(-0.12)
    sigma = sigA*sigB
    # a_T = (aA*aB*aC*aD)*aE*sigma - Split into Terms
    aA = 0.026/((2 * R_t)**0.2)
    aB = ((visc_ns**0.2) * cp_ns)/(praneff_ns**0.6)
    aC = (P_c * g/C_star)**0.8;
    aD = ((2 * R_t)/R_tCurve)**0.1
    aE = ((R_t**2)/(contourR**2))**0.9
    h_g = (aA * aB * aC * aD) * aE * sigma #Convection coefficient
    q_conv = h_g * (T_aw - T_wg) #Heat flux (heat per area)

    # Re-Unit Convert
    h_g = h_g * 144 * 20428.175 #BTU/in^2-s-F -> BTU/ft^2-s-F -> W/m^2K
    q_conv = q_conv * 1634246.235 #BTU/in^2-s -> W/m^2
    T_aw = T_aw / 1.8 #R -> K

    # print(h_g)
    return (h_g, q_conv, T_aw)
