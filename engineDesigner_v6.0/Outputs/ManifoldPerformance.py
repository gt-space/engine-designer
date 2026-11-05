import numpy as np
import matplotlib.pyplot as plt


def plotManifoldResults(contour, opt_results, cad_data):
    fig = plt.figure(figsize=(14, 9))
    fig.suptitle(
        f"Single-Inlet Tapered Manifold Optimization | "
        f"Inlet Split Angle: {opt_results['inlet_angle_deg']:.1f}° | "
        f"Max Maldistribution: {opt_results['maldistribution_pct']:.2f}% | "
        f"Total System dP: {opt_results['dP_total_Pa'] / 6894.75729:.1f} psi",
        fontsize=13, fontweight="bold"
    )

    theta_deg = opt_results["theta_stations_deg"]
    
    # Panel 1: Manifold Tapered Radius & Area Schedule
    ax1 = fig.add_subplot(2, 2, 1)
    r_mm = opt_results["r_manifold_profile"] * 1e3
    t_mm = opt_results["t_wall_profile"] * 1e3
    ctrl_theta = np.linspace(0.0, 180.0, len(opt_results["r_control_points"]))
    ctrl_r_mm = opt_results["r_control_points"] * 1e3
    
    ax1.plot(theta_deg, r_mm, "b-", lw=2.2, label="Inner Flow Radius $r_{in}(\\theta)$")
    ax1.plot(theta_deg, r_mm + t_mm, "k--", lw=1.8, label="Outer Structural Wall $r_{out}(\\theta)$")
    ax1.scatter(ctrl_theta, ctrl_r_mm, color="red", zorder=5, s=50, label="Spline Control Points")
    ax1.set_xlabel("Circumferential Angle $\\theta$ [deg] (0° = Inlet, 180° = End)")
    ax1.set_ylabel("Manifold Section Radius [mm]")
    ax1.set_title("1. Optimized Manifold Taper Profile")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="upper right")

    # Panel 2: Channel Mass Flow Uniformity
    ax2 = fig.add_subplot(2, 2, 2)
    mdot_ch = opt_results["mdot_channels"] * 1e3  # g/s
    mdot_mean = np.mean(mdot_ch)
    ax2.plot(theta_deg, mdot_ch, "g-o", markersize=4, lw=1.8, label="Channel $\\dot{m}_i$")
    ax2.axhline(mdot_mean, color="k", linestyle=":", label=f"Target Mean ({mdot_mean:.2f} g/s)")
    ax2.fill_between(theta_deg, mdot_mean * 0.99, mdot_mean * 1.01, color="green", alpha=0.15, label="±1% Uniformity Band")
    ax2.set_xlabel("Channel Angular Position $\\theta$ [deg]")
    ax2.set_ylabel("Channel Mass Flow Rate [g/s]")
    ax2.set_title("2. Cooling Channel Flow Distribution")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="best")

    # Panel 3: Static Pressure Balance Around the Ring
    ax3 = fig.add_subplot(2, 2, 3)
    P_stat_psi = opt_results["P_static_manifold"] / 6894.75729
    ax3.plot(theta_deg, P_stat_psi, "m-s", markersize=4, lw=1.8)
    ax3.set_xlabel("Circumferential Angle $\\theta$ [deg]")
    ax3.set_ylabel("Manifold Static Pressure [psi]")
    ax3.set_title("3. Circumferential Static Pressure (Friction vs Regain)")
    ax3.grid(True, alpha=0.3)

    # Panel 4: 3D Engine Contour + Tapered Manifold Ring Preview
    ax4 = fig.add_subplot(2, 2, 4, projection="3d")
    X_in, Y_in, Z_in = cad_data["mesh_inner"]
    ax4.plot_surface(X_in * 1e3, Y_in * 1e3, Z_in * 1e3, cmap="coolwarm", alpha=0.85, edgecolor="none")
    
    z_nozzle = contour["z"][-40:] * 1e3
    r_nozzle = contour["r"][-40:] * 1e3
    phi_noz = np.linspace(0, 2 * np.pi, 40)
    Z_N, PHI_N = np.meshgrid(z_nozzle, phi_noz)
    R_N = np.tile(r_nozzle, (len(phi_noz), 1))
    ax4.plot_wireframe(R_N * np.cos(PHI_N), R_N * np.sin(PHI_N), Z_N, color="gray", alpha=0.3, lw=0.6)
    
    ax4.set_xlabel("X [mm]")
    ax4.set_ylabel("Y [mm]")
    ax4.set_zlabel("Axial Z [mm]")
    ax4.set_title("4. 3D Tapered Manifold Loft on Nozzle Exit")

    plt.tight_layout(rect=[0, 0, 1, 0.95])