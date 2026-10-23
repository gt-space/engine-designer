import matplotlib.pyplot as plt

def plotGeometry(contour, channel):

    z_in = contour["z"] * 39.3701
    r_in = contour["r"] * 39.3701

    N_ch = channel["N_ch"]

    w_ch = channel["w_ch"] * 1000 # m -> mm
    h_ch = channel["h_ch"] * 1000  # m -> mm
    w_fin = channel["w_fin"] * 1000 # m -> mm
    tw = channel["tw"] * 1000  # m -> mm

    AR = channel["AR_ch"]
    A_ch = channel["A_ch"] * 1e6 # m^2 -> mm^2


    # -----------------------------
    # Create figure with 4 subplots
    # -----------------------------
    fig, axs = plt.subplots(2, 2, figsize=(14, 9))

    fig.suptitle("Engine and Cooling Channel Geometry")


    # -----------------------------
    # 1. Engine Contour
    # -----------------------------
    axs[0, 0].plot(z_in, r_in)
    axs[0, 0].set_xlabel("Axial Position, z [in]")
    axs[0, 0].set_ylabel("Radius, r [in]")
    axs[0, 0].set_title("Engine Contour")
    axs[0, 0].grid()
    axs[0, 0].axis("equal")


    # -----------------------------
    # 2. Channel Dimensions
    # -----------------------------
    axs[0, 1].plot(z_in, w_ch, label="Channel Width")
    axs[0, 1].plot(z_in, h_ch, label="Channel Height")
    axs[0, 1].plot(z_in, w_fin, label="Fin Width")
    axs[0, 1].plot(z_in, tw, label="Wall Thickness")

    axs[0, 1].set_xlabel("Axial Position, z [in]")
    axs[0, 1].set_ylabel("Dimension [mm]")
    axs[0, 1].set_title(f"Cooling Channel Geometry ({N_ch} Channels)")
    axs[0, 1].set_ylim(bottom=0)
    axs[0, 1].legend()
    axs[0, 1].grid()


    # -----------------------------
    # 3. Aspect Ratio
    # -----------------------------
    axs[1, 0].plot(z_in, AR)

    axs[1, 0].set_xlabel("Axial Position, z [in]")
    axs[1, 0].set_ylabel("Aspect Ratio, h/w")
    axs[1, 0].set_title("Cooling Channel Aspect Ratio")
    axs[1, 0].grid()


    # -----------------------------
    # 4. Channel Area
    # -----------------------------
    axs[1, 1].plot(z_in, A_ch)

    axs[1, 1].set_xlabel("Axial Position, z [in]")
    axs[1, 1].set_ylabel("Channel Area [mm²]")
    axs[1, 1].set_title("Cooling Channel Cross-Sectional Area")
    axs[1, 1].grid()


    # Keep everything from overlapping
    plt.tight_layout(rect=[0, 0, 1, 0.96])