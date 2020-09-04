def getProps(T_cb):
    #jetAProps Solve for fluid properties of Jet-A given temperature and
    #pressure - No correlations are made to pressure but it will
    #be included for future development as compressibility effects are
    #considered

    # Properties derived from 2 pt. linear correlation using DTIC A132106
    # "Aviation Fuel Properties"

    #Solve for following properties:
    # - Density
    # - Viscosity
    # - Conductivity
    # - Specific heat

    # https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=904848 See page 73
    Tc_cel = T_cb - 273.15 #Coolant temperature in celcius.

    # Density
    # Temperature in Celcius
    # rho_90C = 750 kg/m^3, rho_-40C = 848 kg/m^3, Pg. 22
    rhoSlope = -0.753846
    rho_c = rhoSlope * (Tc_cel + 40) + 848 #kg/m^3

    # Specific Heat
    # Temperature in Celcius
    # Cp_40C = 2.04 kj/kgK, Cp_160C = 2.56 kj/kgK, Pg. 55
    CpSlope = 0.004333333
    C_pc = 1000 * (CpSlope * (Tc_cel - 40) + 2.04) #J/kgK

    # Conductivity
    # Temperature in Celcius
    # Cond_20C = 0.1150 W/mK, Cond_100C = 0.1010 W/mK, Pg. 57
    condSlope = -0.000175
    cond_c = condSlope * (Tc_cel - 20) + 1.1150 #W/mK

    viscK_c = 2.3498 / 1000000 #Specific val for ex. Will need to set up interpolation system (mm^2/s => m^2/s)

    return (rho_c, C_pc, cond_c, viscK_c)
