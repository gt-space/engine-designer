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
    Tc_cel = T_cb - 273.15 #Coolant temperature in celcius.

    # Boiling point at 1 atm (conservative since BP increases with pressure)
    T_boil = 462.039 # (K)

    # Density
    # Temperature in Celcius
    # rho_288.7K = 0.811 kg/L
    rhoSlope = -0.753846 # From old data. Would be nice to update
    rho_c = rhoSlope * (T_cb - 288.7) + 848 #kg/m^3

    # Specific Heat
    # Temperature in Kelvin
    # Cp_250K = 1.73 J/gK, Cp_450K = 2.55 J/gK, Pg. 85
    CpSlope = (2.55-1.73)/(450-250)
    C_pc = 1000 * (CpSlope * (T_cb - 250) + 1.73) #J/gK => J/kgK

    # Cp_250K = 1.73 kj/kgK, Cp_450K = 2.55 kj.kgK
    CpSlope = 0.0041

    # Conductivity
    # Temperature in Celcius
    # Cond_20C = 0.1150 W/mK, Cond_100C = 0.1010 W/mK, Pg. 82
    condSlope = -0.000175
    cond_c = condSlope * (Tc_cel - 20) + 1.1150 #W/mK

    # Kinematic viscosity
    viscK_c = (11.5 / 1000000) * math.exp(-(T_cb-238.706)/60) # Taken at -30 F system (mm^2/s => m^2/s)
        # 60 comes from decay rate of this data: # https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=904848 (page 73)

    return (rho_c, C_pc, cond_c, viscK_c)
