import numpy as np
from Functions.CoolantProperties import (get_coolant_h, get_coolant_props)
from Functions.HeatTransfer import (coolant_htc, bartz_hg, fin_efficiency)
from Functions.PressureDrop import dPdarcyWeisbach, dPfluidAcceleration
from EngineParameters import mdot, MR, P_cool_out_desired
from Geometry.Contour import getCurrentEngineParameters
from Solvers.HotGasFlowSolver import getGasPropertiesCEA



def CoolantChannelSolver(T_in_coolant, contour, channel, hot_gas):
    gas = getGasPropertiesCEA()
    geometry_params = getCurrentEngineParameters()

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

    # Helper to solve radial 1D thermal resistance network at any station i
    def solve_station_thermal(i, Tc, hc, ds_local, Twg_init, Twc_init):
        Ag = (w_ch[i] + w_fin[i]) * ds_local
        Twg_old, Twc_old = Twg_init, Twc_init

        for _ in range(400):
            Tw_avg = 0.5 * (Twg_old + Twc_old)
            kw = material.get("thermal_conductivity", T=Tw_avg)
            hg = bartz_hg(A[i], M[i], Twg_old)

            N_fin = fin_efficiency(hc, kw, h_ch[i], w_fin[i])
            A_eff = (w_ch[i] + 2.0 * N_fin * h_ch[i]) * ds_local

            R_g = 1.0 / (hg * Ag)
            R_wall = tw[i] / (Ag * kw)
            R_c = 1.0 / (hc * A_eff)
            Q = (Taw[i] - Tc) / (R_g + R_wall + R_c)

            Twg_new = Taw[i] - Q / (hg * Ag)
            Twc_new = Tc + Q / (hc * A_eff)

            if abs(Twg_new - Twg_old) < 0.01 and abs(Twc_new - Twc_old) < 0.01:
                return Q, Q / Ag, Twg_new, Twc_new

            Twg_old, Twc_old = Twg_new, Twc_new

        return Q, Q / Ag, Twg_old, Twc_old

    #make a guess for inlet pressure
    P_in_guess = P_cool_out_desired + 40.0 * 6894.75729 #guessing 400psi + convert to Pa

    for _ in range(400):
        P_cool_array[N - 1] = P_in_guess
        h_cool_array[N - 1] = get_coolant_h(T_in_coolant, P_in_guess)

        for i in range(N - 1, 0, -1):
            coolant_props = get_coolant_props(P_cool_array[i], h_cool_array[i])
            Tc = coolant_props["T"]
            rho_up = coolant_props["rho"]
            mu_c = coolant_props["mu"]
            k_c = coolant_props["k"]
            Pr_c = coolant_props["Pr"]

            T_cool_array[i] = Tc
            Vc = mdot_ch / (rho_up * A_ch[i])
            Re = rho_up * Vc * Dh[i] / mu_c
            hc, _ = coolant_htc(Re, Pr_c, k_c, roughness, Dh[i])

            Twg_init = 800.0 if i == N - 1 else Twg_array[i + 1]
            Twc_init = 700.0 if i == N - 1 else Twc_array[i + 1]

            Q, q_flux, Twg_conv, Twc_conv = solve_station_thermal(
                i, Tc, hc, ds[i - 1], Twg_init, Twc_init
            )
            q_flux_array[i] = q_flux
            Twg_array[i] = Twg_conv
            Twc_array[i] = Twc_conv

            # Update downstream enthalpy (station i-1)
            h_cool_array[i - 1] = h_cool_array[i] + Q / mdot_ch

            # Evaluate downstream density rho_down for true thermal acceleration dP
            props_down = get_coolant_props(P_cool_array[i], h_cool_array[i - 1])
            rho_down = props_down["rho"]

            dP_friction = dPdarcyWeisbach(Re, rho_up, Vc, ds[i - 1], roughness, Dh[i])
            dP_acceleration = dPfluidAcceleration(A_ch[i - 1], A_ch[i], rho_down, rho_up, mdot_ch)

            P_cool_array[i - 1] = P_cool_array[i] - dP_friction - dP_acceleration
            if P_cool_array[i - 1] <= 0:
                raise ValueError("Coolant pressure became non-positive")

        # Populate station i = 0 (Injector exit) so arrays have no zeros at index 0
        props_0 = get_coolant_props(P_cool_array[0], h_cool_array[0])
        T_cool_array[0] = props_0["T"]
        Vc_0 = mdot_ch / (props_0["rho"] * A_ch[0])
        Re_0 = props_0["rho"] * Vc_0 * Dh[0] / props_0["mu"]
        hc_0, _ = coolant_htc(Re_0, props_0["Pr"], props_0["k"], roughness, Dh[0])
        _, q_flux_array[0], Twg_array[0], Twc_array[0] = solve_station_thermal(
            0, props_0["T"], hc_0, ds[0], Twg_array[1], Twc_array[1]
        )

        # Check outlet pressure convergence
        P_cool_out = P_cool_array[0]
        error_Pressure = P_cool_out_desired - P_cool_out
        if abs(error_Pressure) < 100.0:
            break

        P_in_guess += error_Pressure

    return {
        "P_cool": P_cool_array,
        "h_cool": h_cool_array,
        "T_cool": T_cool_array,
        "q_flux": q_flux_array,
        "Twg": Twg_array,
        "Twc": Twc_array
    }