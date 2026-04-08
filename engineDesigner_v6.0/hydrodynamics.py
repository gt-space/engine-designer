def solve_couette(Gamm_L, liquid_props, P_D, C_f)
    #calculate fillm thickness(delta_L)
    numerator = 2 * liquid_props.mu_L * Gamm_L
    denominator = P_D * C_f * liquid_props.rho_l
    delta_L = (numerator / denominator)


    #calculate liquid surface velocity (U_ls)
    U_ls = (C_f * P_D * delta_L) / (2 * liquid_props.mu_l)

    return delta_L, U_ls




def calculate_critical_gas_velocity(gas_props, liquid_props, delta_L, lambd, sigma, g_0, delta_rho):
    #calculate critical gas velocity (U_crit)

    #liquid film thickness (h_l) and gas film thickness (h_g)
    h_l = delta_L
    h_g = #chamber diameter - h_l  can be approximated as chamber diameter for thin films
    import math
    #liquid and gas factors (a_l and a_g), where k is the wavenumber of the disturbance
    k = 2 * math.pi / lambd #wavelength of distubance and is taken from wave characterization,
                               #this is peridoxical but trust me bro, this should either be calculated
                               #using some set of initial values in the wave characterization function 
                               #or can be treated as equal to film thickness at injection point
    a_l = math.coth(k*h_l)
    a_g = math.coth(k*h_g)
    
    rho_ratio = liquid_props.rho_l / gas_props.rho_g
    Beta = math.sqrt(sigma/(g_0 * liquid_props.rho_l)) #surface tension factor, where sigma is liquid surface tension

    term1 = (g_0 *(a_l *liquid_props.mu_l + a_g * gas_props.mu_g)^2)/(a_l * a_g^2 * liquid_props.mu_g^2 + rho_ratio *a_g * a_l^2 * liquid_props.mu_l^2)
    
    term2 = Beta^2 * k + (1-rho_ratio)/ k
    U_crit = math.sqrt(term1 * term2)
    
    return U_crit  




def  check_stability_roughness(U_g, U_ls, delta_L, gas_props, liquid_props, G_g, lambd, sigma, g_0, delta_rho, D, T_g, T_L_sat, k_nu, k_c):
    import math
    #1. Check PKHI
    U_crit = calculate_critical_gas_velocity(gas_props, liquid_props, delta_L, lambd, sigma, g_0, delta_rho)

    if U_g > U_crit:
        #waves exist, characterize them
        #wave number (N_mu)
        N_mu = liquid_props.mu_l / (liquid_props.rho_l * sigma * math.sqrt(sigma/(g_0 * delta_rho)))
        #wave constant (C_w)
        if N_mu == 1/15 or N_mu < 1/15:
            C_w = .028 * N_mu**(-4/5)
        elif N_mu > 1/15:
            C_w = .25
        else:
            print('Error: Invalid N_mu value')

        #Eotvos number (E_o), where D is the chamber diameter
        E_o = (g_0 * D**2 *(liquid_props.rho_l - gas_props.rho_g))/sigma 

        #wave frequency (nu), where k_nu is the proportionality constant ( equal to 1 according to bernal)
        nu = (math.sqrt(U_g * U_ls) / D) * k_nu * gas_props.Re_g**(.53) * liquid_props.Re_l**(-.48) * E_o**(.27) * (liquid_props.rho_l / gas_props.rho_g)**(.14) * C_w**(.68)

        #wave celerity (c), where k_c is the proportionality constant (equal to 1 according to bernal)
        term1 = (math.sqrt(gas_props.rho_g) * U_g + math.sqrt(liquid_props.rho_l) * U_ls) / (math.sqrt(gas_props.rho_g) + math.sqrt(liquid_props.rho_l))
        term2 = k_c * gas_props.Re_g**(-.38) * liquid_props.Re_l**(.16) * C_w**(-.13)
        c = term1 * term2

        lambd = c / nu

        #mass flux of gas (G_g) and Surface Roughness Factor (Fr)
        G_g = gas_props.rho_g * U_g
        Fr = 14.1 * gas_props.rho_g ** (.4) / (G_g ** (.8) * (T_g / T_L_sat) ** (.2))

        return c, nu, lambd, Fr
    else:
        print('No waves exist, flow is stable')
        return None, None, None, None   
    
