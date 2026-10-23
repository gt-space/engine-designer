import matplotlib.pyplot as plt

def plotHotGasFlow(contour, hot_gas):

    z = contour["z"] * 39.3701
    M = hot_gas["M"]
    P = hot_gas["P"] / 6894.75729 #psi
    T = hot_gas["T"]
    Taw = hot_gas["Taw"]
    V = hot_gas["V"]
    rho = hot_gas["rho"]

    fig, axs = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle("Hot Gas Flow Properties")

    # Row 1, Col 1
    axs[0, 0].plot(z, M)
    axs[0, 0].set_title("Mach Number")
    axs[0, 0].set_xlabel("Axial Position [in]")
    axs[0, 0].set_ylabel("M")
    axs[0, 0].grid()

    # Row 1, Col 2
    axs[0, 1].plot(z, P)
    axs[0, 1].set_title("Static Pressure")
    axs[0, 1].set_xlabel("Axial Position [in]")
    axs[0, 1].set_ylabel("Pressure [Psi]")
    axs[0, 1].grid()

    # Row 1, Col 3
    axs[0, 2].plot(z, T)
    axs[0, 2].set_title("Static Temperature")
    axs[0, 2].set_xlabel("Axial Position [in]")
    axs[0, 2].set_ylabel("Temperature [K]")
    axs[0, 2].grid()

    # Row 2, Col 1
    axs[1, 0].plot(z, Taw)
    axs[1, 0].set_title("Adiabatic Wall Temperature")
    axs[1, 0].set_xlabel("Axial Position [in]")
    axs[1, 0].set_ylabel("Temperature [K]")
    axs[1, 0].grid()

    # Row 2, Col 2
    axs[1, 1].plot(z, V)
    axs[1, 1].set_title("Gas Velocity")
    axs[1, 1].set_xlabel("Axial Position [in]")
    axs[1, 1].set_ylabel("Velocity [m/s]")
    axs[1, 1].grid()

    # Row 2, Col 3
    axs[1, 2].plot(z, rho)
    axs[1, 2].set_title("Density")
    axs[1, 2].set_xlabel("Axial Position [in]")
    axs[1, 2].set_ylabel("Density [kg/m³]")
    axs[1, 2].grid()

    plt.tight_layout()