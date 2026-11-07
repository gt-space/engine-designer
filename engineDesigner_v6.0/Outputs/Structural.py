import numpy as np
import matplotlib.pyplot as plt


def plotStructuralSizing(contour, struct_results):
    z_in = contour["z"] * 39.3701
    r_inner_mm = contour["r"] * 1e3

    liner = struct_results["liner"]
    shell = struct_results["shell"]
    ch = struct_results["channel_updated"]

    fig, axs = plt.subplots(2, 2, figsize=(14, 9))
    fig.suptitle(
        f"Regen Liner ({liner['material_name']}) & Outer Shell ({shell['shell_material_name']}) Sizing | "
        f"Chamber Mass: {struct_results['mass_total_chamber_kg']:.2f} kg",
        fontsize=13, fontweight="bold"
    )

    # 1. Inner Liner Thickness Comparison
    axs[0, 0].plot(z_in, liner["t_liner_structural"] * 1e3, "b--", lw=1.8, label="Structural Thickness $t_{w,struct}$")
    axs[0, 0].plot(z_in, liner["t_liner_thermal"] * 1e3, "r:", lw=1.8, label="Thermal/Ratcheting $t_{w,thermal}$")
    axs[0, 0].plot(z_in, liner["t_liner_selected"] * 1e3, "k-", lw=2.5, label="Selected Priority $\\max(t_{struct}, t_{th})$")
    axs[0, 0].set_title(f"1. Inner Liner Thickness ({liner['material_name']})")
    axs[0, 0].set_xlabel("Axial Position, z [in]")
    axs[0, 0].set_ylabel("Liner Thickness [mm]")
    axs[0, 0].legend(loc="best")
    axs[0, 0].grid(True, alpha=0.3)

    # 2. Liner Selected Thickness vs. 1D Fourier Conduction Ceiling
    axs[0, 1].plot(z_in, liner["t_liner_selected"] * 1e3, "k-", lw=2.2, label="Selected Liner $t_w$")
    axs[0, 1].plot(z_in, np.minimum(liner["t_liner_max_fourier"] * 1e3, 15.0), "r--", lw=1.8,
                   label="Max Fourier Thermal Ceiling $t_{w,max}$")
    axs[0, 1].fill_between(
        z_in, liner["t_liner_selected"] * 1e3,
        np.minimum(liner["t_liner_max_fourier"] * 1e3, 15.0),
        color="green", alpha=0.12, label="Safe Conduction Margin"
    )
    axs[0, 1].set_title("2. Liner Conduction Overheat Margin")
    axs[0, 1].set_xlabel("Axial Position, z [in]")
    axs[0, 1].set_ylabel("Thickness [mm]")
    axs[0, 1].legend(loc="best")
    axs[0, 1].grid(True, alpha=0.3)

    # 3. Outer Closeout Shell Thickness Comparison
    axs[1, 0].plot(z_in, shell["t_shell_structural"] * 1e3, "b--", lw=1.8, label="Structural Vessel $t_{shell,struct}$")
    axs[1, 0].plot(z_in, shell["t_shell_thermal"] * 1e3, "r:", lw=1.8, label="Thermal Expansion $t_{shell,thermal}$")
    axs[1, 0].plot(z_in, shell["t_shell_selected"] * 1e3, "k-", lw=2.5, label="Selected Priority $\\max(t_{struct}, t_{th})$")
    axs[1, 0].set_title(f"3. Outer Shell Thickness After Channels ({shell['shell_material_name']})")
    axs[1, 0].set_xlabel("Axial Position, z [in]")
    axs[1, 0].set_ylabel("Outer Shell Thickness [mm]")
    axs[1, 0].legend(loc="best")
    axs[1, 0].grid(True, alpha=0.3)

    # 4. Complete Radial CAD Wall Stack-Up
    r_liner_out_mm = r_inner_mm + ch["tw"] * 1e3
    r_ch_out_mm = shell["r_shell_inner"] * 1e3
    r_shell_out_mm = shell["r_shell_outer"] * 1e3

    axs[1, 1].plot(z_in, r_inner_mm, "r-", lw=1.8, label="Hot-Gas Wall $r_{gas}$")
    axs[1, 1].plot(z_in, r_liner_out_mm, "orange", lw=1.5, label="Liner Outer / Channel Floor")
    axs[1, 1].plot(z_in, r_ch_out_mm, "c--", lw=1.5, label="Channel Roof / Shell Inner $r_{shell,in}$")
    axs[1, 1].plot(z_in, r_shell_out_mm, "k-", lw=2.2, label="Outer Structural Shell $r_{shell,out}$")
    axs[1, 1].fill_between(z_in, r_inner_mm, r_liner_out_mm, color="orange", alpha=0.35)
    axs[1, 1].fill_between(z_in, r_liner_out_mm, r_ch_out_mm, color="cyan", alpha=0.25)
    axs[1, 1].fill_between(z_in, r_ch_out_mm, r_shell_out_mm, color="gray", alpha=0.45)

    axs[1, 1].set_title("4. Full 3D CAD Radial Wall Stack-Up")
    axs[1, 1].set_xlabel("Axial Position, z [in]")
    axs[1, 1].set_ylabel("Radius [mm]")
    axs[1, 1].legend(loc="best")
    axs[1, 1].grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.95])