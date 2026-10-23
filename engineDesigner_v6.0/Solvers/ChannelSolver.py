import numpy as np
from Functions.CoolantProperties import (get_coolant_h, get_coolant_props)
from Functions.HeatTransfer import (coolant_htc, bartz_hg, fin_efficiency)
from Functions.PressureDrop import dPdarcyWeisbach, dPfluidAcceleration
from EngineParameters import mdot, MR, P_cool_out_desired



def CoolantChannelSolver(T_in_coolant, contour, channel, hot_gas):

    #Get contour geometry
    N = contour["N"]
    A = contour["A"]
    ds = contour["ds"]

    #Get channel geometry
    N_ch = channel["N_ch"]
    w_ch = channel["w_ch"]
    w_fin = channel["w_fin"]
    h_ch = channel["h_ch"]
    A_ch = channel["A_ch"]
    Dh = channel["Dh"]
    tw = channel["tw"]
    roughness = channel["roughness"]
    material = channel["material"]

    #get hot gas flow
    M = hot_gas["M"]
    Taw = hot_gas["Taw"]

    #get channel flow rate
    mdot_f = mdot / (1 + MR)
    mdot_ch = mdot_f/N_ch

    #Setup Arrays to calculate thermo-fluid properties at every station along the channel
    P_cool_array = np.zeros(N) #pressure along channel
    h_cool_array = np.zeros(N) #enthalpy along channel
    T_cool_array = np.zeros(N) #bulk temperature along channel 
    q_flux_array = np.zeros(N) #heat flux along channel 
    Twg_array = np.zeros(N) #hot-gas side wall temperature along channel
    Twc_array = np.zeros(N) # cold-gas side wall temperature along channel

    #make a guess for inlet pressure
    P_in_guess = 400*6895 #guessing 400psi + convert to Pa

    for pressure_iteration in range(400):
        P_cool_array[N-1] = P_in_guess
        h_cool_array[N-1] = get_coolant_h(T_in_coolant, P_in_guess)

        #range N is the number of axial points on the engine
        for i in range(N-1,0,-1):

            #get properties of the fluid based on pressure and enthalpy of the fluid
            coolant_props = get_coolant_props(P_cool_array[i], h_cool_array[i])
            Tc, rho_c, mu_c, k_c, cp_c, Pr_c = coolant_props.values()

            T_cool_array[i] = Tc
            Vc = mdot_ch/(rho_c*A_ch[i])
            Re = rho_c*Vc*Dh[i]/mu_c
            hc, Nu = coolant_htc(Re, Pr_c, k_c, roughness, Dh[i])

            Ag = (w_ch[i]+w_fin[i])*ds[i-1] #hot gas wall area
        

            #Make initial guess for temperature
            if i == N-1: #use these values for first guess 
                Twg_old = 800.0
                Twc_old = 700.0
            else: # Use previously solved station as initial guess
                Twg_old = Twg_array[i+1]
                Twc_old = Twc_array[i+1]
            
            for iteration in range(400):
                Tw_avg = (Twg_old+Twc_old)/2 #using average wall temp for k wall
                kw = material.get("thermal_conductivity", Tw_avg) #get thermal conductivity based on average wall temp
                hg = bartz_hg(A[i], M[i], Twg_old) #calculate hg using bartz

                N_fin = fin_efficiency(hc, kw, h_ch[i], w_fin[i]) #calculate fin efficiency
                A_eff = (w_ch[i]+2*N_fin*h_ch[i])*ds[i-1]


                R_g = 1/(hg*Ag); #resistance, hot gas
                R_wall = tw[i]/(Ag*kw) #resistance, wall
                R_c = 1/(hc*A_eff) #resistance, coolant (with fin efficiency)
                Q = (Taw[i]-Tc)/(R_g+R_wall+R_c) #Calculate heat transfer rate
                Twg_new = Taw[i]-Q/(hg*Ag)
                Twc_new = Tc+Q/(hc*A_eff)

                #Check for convergence
                diff_Twg = abs(Twg_new-Twg_old)
                diff_Twc = abs(Twc_new-Twc_old)

                #Set new guess to calculated wall temperatures 
                Twg_old = Twg_new
                Twc_old = Twc_new

                #break out of the iteration if converged
                if diff_Twg <0.01 and diff_Twc <0.01:
                    q_flux_array[i] = Q/Ag
                    Twg_array[i] = Twg_new
                    Twc_array[i] = Twc_new
                    break

            #update enthalpy and pressure based on inner loop iteration values
            #calculate pressure error
            h_cool_array[i-1] = (h_cool_array[i] + Q/mdot_ch)
            dP_friction = dPdarcyWeisbach(Re, rho_c, Vc, ds[i-1], roughness, Dh[i])
            dP_acceleration = dPfluidAcceleration(A_ch[i-1], A_ch[i], rho_c, mdot_ch)
            #dP_areaChange = dPareaChange()
            #dP_fluidAcceleration = dPfluidAcceleration()

            P_cool_array[i-1] = P_cool_array[i] - dP_friction - dP_acceleration
            if P_cool_array[i] <= 0:
                raise ValueError("Coolant pressure became non-positive")

        #Check for convergence
        P_cool_out = P_cool_array[0] #calculated outlet pressure

        #compare calculated outlet pressure to actual outlet pressure
        error_Pressure = P_cool_out_desired - P_cool_out

        #break out of the iteration if true
        if abs(error_Pressure) <100: #Pa
            break

        #update inlet pressure guess if not converged
        P_in_guess = P_in_guess + error_Pressure

    return {
    "P_cool": P_cool_array,
    "h_cool": h_cool_array,
    "T_cool": T_cool_array,
    "q_flux": q_flux_array,
    "Twg": Twg_array,
    "Twc": Twc_array 
    }