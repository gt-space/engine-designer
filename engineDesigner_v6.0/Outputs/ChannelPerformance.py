import matplotlib.pyplot as plt

def plotChannelPerformance(contour, hot_gas, coolant_flow):

    # Axial position
    z_in = contour["z"] * 39.3701

    # Hot gas
    Taw = hot_gas["Taw"]

    # Coolant / wall solution
    Twg = coolant_flow["Twg"]
    Twc = coolant_flow["Twc"]
    T_cool = coolant_flow["T_cool"]
    P_cool = coolant_flow["P_cool"] / 6894.75729   # Pa -> psi
    q_flux = coolant_flow["q_flux"] / 1e6          # W/m² -> W/cm²
    h_cool = coolant_flow["h_cool"] / 1000         # J/kg -> kJ/kg


    # Create figure
    fig, axs = plt.subplots(2, 2, figsize=(14, 9))

    fig.suptitle("Cooling Channel Performance")


    # ------------------------------------------------
    # Temperatures
    # ------------------------------------------------
    axs[0, 0].plot(
        z_in[1:], Taw[1:],
        label=r"Adiabatic Wall Temperature $T_{aw}$")

    axs[0, 0].plot(
        z_in[1:], Twg[1:],
        label=r"Gas-Side Wall Temperature $T_{wg}$")

    axs[0, 0].plot(
        z_in[1:], Twc[1:],
        label=r"Coolant-Side Wall Temperature $T_{wc}$")

    axs[0, 0].plot(
        z_in[1:], T_cool[1:],
        label=r"Bulk Coolant Temperature $T_c$")

    axs[0, 0].set_title("Temperatures Along Engine")
    axs[0, 0].set_xlabel("Axial Position [in]")
    axs[0, 0].set_ylabel("Temperature [K]")
    axs[0, 0].legend()
    axs[0, 0].grid()


    # ------------------------------------------------
    # Coolant pressure
    # ------------------------------------------------
    axs[0, 1].plot(z_in[1:], P_cool[1:])

    axs[0, 1].set_title("Coolant Pressure Along Channel")
    axs[0, 1].set_xlabel("Axial Position [in]")
    axs[0, 1].set_ylabel("Coolant Pressure [psi]")
    axs[0, 1].grid()


    # ------------------------------------------------
    # Heat flux
    # ------------------------------------------------
    axs[1, 0].plot(z_in[1:], q_flux[1:])

    axs[1, 0].set_title("Gas-Side Heat Flux")
    axs[1, 0].set_xlabel("Axial Position [in]")
    axs[1, 0].set_ylabel("Heat Flux [MW/m²]")
    axs[1, 0].grid()


    # ------------------------------------------------
    # Coolant enthalpy
    # ------------------------------------------------
    axs[1, 1].plot(z_in[1:], h_cool[1:])

    axs[1, 1].set_title("Coolant Specific Enthalpy")
    axs[1, 1].set_xlabel("Axial Position [in]")
    axs[1, 1].set_ylabel("Specific Enthalpy [kJ/kg]")
    axs[1, 1].grid()


    plt.tight_layout(rect=[0, 0, 1, 0.95])