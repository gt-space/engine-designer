
def Bartz(engineProps, bartzData, i):
    # unpack data
    T_wg = bartzData[0]
    T_c0 = bartzData[1]
    R_t = bartzData[2]
    P_c = bartzData[3]
    C_star = bartzData[4]
    R_tCurve = bartzData[5]
    praneff_ns = bartzData[6]
    visc_ns = bartzData[7]
    cp_ns = bartzData[8]
    g = bartzData[9]

    # Run Bartz Correlation
    gam = engineProps[i, 15] #Index gamma; No conversion needed
    mach = engineProps[i, 5]
    contourR = engineProps[i, 0] * 39.3701 #Convert to Inches
    T_aw = T_c0 * (1 + (praneff_ns**(1/3)) * ((gam - 1)/2) * (mach**2))/(1 + ((gam - 1)/2) * (mach**2))
    #Sigma = sigA * sigB - Split into Terms for readability
    sigA = (0.5 * (T_wg/T_c0) * (1 + ((gam - 1)/2) * mach**2) + 0.5)**(-0.68)
    sigB = (1 + ((gam - 1)/2) * (mach**2))**(-0.12);
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

    return (h_g, q_conv, T_aw)
