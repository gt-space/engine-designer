def getProps(T_cb):
    #Solve for fluid properties of Jet-A given temperature and
    #pressure - No correlations are made to pressure but it will
    #be included for future development as compressibility effects are
    #considered

    #Solve for following properties:
    # - Density
    # - Viscosity
    # - Conductivity
    # - Specific heat

    import math

    # https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=904848 See page 73

    # Boiling point at 1 atm (conservative since BP increases with pressure)
    T_boil = 462.039 # (K)

    # Density
    # Temperature in Celcius
    # rho_270K = 0.835 kg/L
    # rho_470K = 0.685 kg/L
    rhoSlope = (0.685-0.835)/(470-270) # From old data. Would be nice to update
    rho_c = (rhoSlope * (T_cb - 270) + 0.835) * 1000 #kg/m^3

    # Specific Heat
    # Temperature in Kelvin
    # Cp_250K = 1.73 J/gK, Cp_450K = 2.55 J/gK, Pg. 85
    CpSlope = (2.55-1.73)/(450-250)
    C_pc = 1000 * (CpSlope * (T_cb - 250) + 1.73) #J/gK => J/kgK

    # Cp_250K = 1.73 kj/kgK, Cp_450K = 2.55 kj.kgK
    CpSlope = 0.0041

    # Conductivity
    # Temperature in Celcius
    # Cond_300K = 0.1150 W/mK, Cond_550K = 0.076 W/mK, Pg. 82
    condSlope = (0.076-0.1150)/(550-300)

    cond_c = condSlope * (T_cb - 300) + 0.1150 #W/mK

    # Kinematic viscosity
    viscK_c =  1.7825 * (1700 * math.exp(-0.026*(T_cb)) + 0.27) / 1000000 # Done with eyeballed curve fit. Should be updated to use a more rigorous method

    return (rho_c, C_pc, cond_c, viscK_c)
