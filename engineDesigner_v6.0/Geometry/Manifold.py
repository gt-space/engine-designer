import numpy as np
from scipy.interpolate import PchipInterpolator

import EngineParameters as ep
from Functions.CoolantProperties import get_coolant_h, get_coolant_props
import Functions.Structural as st
import Optimization.ManifoldOptimization as opt

#Get the geometry data from cooling channels and Contour
def extractManifoldChannelInterface(contour, channel, coolant_flow=None):
    N_ch = channel["N_ch"]
    L_channel = float(np.sum(contour["ds"]))
    A_ch_avg = float(np.mean(channel["A_ch"]))
    Dh_ch_avg = float(np.mean(channel["Dh"]))
    
    # Manifold sits at the nozzle exit (index -1 of contour and channel arrays)
    R_manifold = float(contour["r"][-1] + channel["tw"][-1] + channel["h_ch"][-1])
    z_manifold = float(contour["z"][-1])
    A_ch_inlet = float(channel["A_ch"][-1])
    
    # If ChannelSolver was already run, calibrate f_channel to match its exact dP!
    if coolant_flow is not None:
        dP_ch_actual = float(coolant_flow["P_cool"][-1] - coolant_flow["P_cool"][0])
        mdot_f = ep.mdot / (1.0 + ep.MR)
        mdot_ch_nom = mdot_f / N_ch
        P_in_ch = float(coolant_flow["P_cool"][-1])
        h_in = get_coolant_h(ep.T_cool_in, P_in_ch)
        rho_in = get_coolant_props(P_in_ch, h_in)["rho"]
        V_avg = mdot_ch_nom / (rho_in * A_ch_avg)
        f_equiv = dP_ch_actual / ((L_channel / Dh_ch_avg) * 0.5 * rho_in * V_avg**2)
    else:
        f_equiv = 0.022

    return {
        "N_channels": int(N_ch),
        "R_manifold": R_manifold,
        "z_manifold": z_manifold,
        "A_ch_inlet": A_ch_inlet,
        "A_channel_avg": A_ch_avg,
        "Dh_channel": Dh_ch_avg,
        "L_channel": L_channel,
        "f_channel": float(np.clip(f_equiv, 0.005, 0.25))
    }

#Sets Manifold Limits based on contour and structural properties
def setManifoldSizeLimitsFromContour(contour, channel=None, max_outer_radius_ratio=1.35, min_mfg_radius=0.0025):
    r_array = contour["r"]
    r_max_contour = float(np.max(r_array))
    r_at_manifold = float(r_array[-1])
    
    max_outer_envelope_radius = r_max_contour * max_outer_radius_ratio
    max_radial_protrusion = max(max_outer_envelope_radius - r_at_manifold, 0.015)
    
    return {
        "yield_strength": 450e6,
        "rho_material": 8400.0,
        "safety_factor": 2.0,
        "t_min_mfg": 0.0015,
        "r_min_mfg": min_mfg_radius,
        "max_radial_envelope": max_radial_protrusion
    }

#Load fuel coolant mdot and fluid props
def loadEngineCoolantParameters(coolant_flow=None, mdot_override=None):
    mdot_fuel = mdot_override if mdot_override is not None else (ep.mdot / (1.0 + ep.MR))
    
    d_feed_pipe = getattr(ep, "d_feed_pipe", 0.75 * 0.0254)
    A_feed_pipe = 0.25 * np.pi * d_feed_pipe**2
    
    if coolant_flow is not None:
        P_ch_inlet = float(coolant_flow["P_cool"][-1])
        P_inlet_guess = P_ch_inlet * 1.05
    else:
        P_inlet_guess = ep.P_cool_out_desired * 1.35
        
    h_in = get_coolant_h(ep.T_cool_in, P_inlet_guess)
    props = get_coolant_props(P_inlet_guess, h_in)
    
    return {
        "mdot_total": float(mdot_fuel),
        "rho": float(props["rho"]),
        "mu": float(props["mu"]),
        "roughness": float(ep.channel_roughness),
        "A_feed_pipe": float(A_feed_pipe),
        "P_inlet_guess": float(P_inlet_guess),
        "P_collector_exit": float(ep.P_cool_out_desired)
    }

#Returns numbers for CAD geometry
def build3DCADManifoldLoft(opt_results, channel_geom, n_cad_angles=361, n_section_pts=64):
    theta_opt_rad = opt_results["theta_stations_rad"]
    r_opt = opt_results["r_manifold_profile"]
    t_opt = opt_results["t_wall_profile"]
    
    R_base = channel_geom["R_manifold"]
    z_base = channel_geom["z_manifold"]
    
    r_spline = PchipInterpolator(theta_opt_rad, r_opt)
    t_spline = PchipInterpolator(theta_opt_rad, t_opt)
    
    theta_full = np.linspace(-np.pi, np.pi, n_cad_angles)
    abs_theta = np.abs(theta_full)
    
    r_inner_full = r_spline(abs_theta)
    t_wall_full = t_spline(abs_theta)
    r_outer_full = r_inner_full + t_wall_full
    
    R_centerline = R_base + r_inner_full
    X_center = R_centerline * np.cos(theta_full)
    Y_center = R_centerline * np.sin(theta_full)
    Z_center = np.full_like(theta_full, z_base)
    
    phi = np.linspace(0.0, 2.0 * np.pi, n_section_pts)
    THETA, PHI = np.meshgrid(theta_full, phi, indexing="ij")
    
    R_IN = r_inner_full[:, np.newaxis]
    R_OUT = r_outer_full[:, np.newaxis]
    R_CEN = R_centerline[:, np.newaxis]
    
    X_inner = (R_CEN + R_IN * np.cos(PHI)) * np.cos(THETA)
    Y_inner = (R_CEN + R_IN * np.cos(PHI)) * np.sin(THETA)
    Z_inner = z_base + R_IN * np.sin(PHI)
    
    X_outer = (R_CEN + R_OUT * np.cos(PHI)) * np.cos(THETA)
    Y_outer = (R_CEN + R_OUT * np.cos(PHI)) * np.sin(THETA)
    Z_outer = z_base + R_OUT * np.sin(PHI)
    
    return {
        "theta_full_deg": np.degrees(theta_full),
        "theta_full_rad": theta_full,
        "r_inner_full_m": r_inner_full,
        "r_outer_full_m": r_outer_full,
        "t_wall_full_m": t_wall_full,
        "R_centerline_m": R_centerline,
        "centerline_xyz": np.column_stack((X_center, Y_center, Z_center)),
        "mesh_inner": (X_inner, Y_inner, Z_inner),
        "mesh_outer": (X_outer, Y_outer, Z_outer)
    }

#Export to CSV
def exportManifoldGeometryToCSV(cad_data, opt_results, filename="Geometry/Manifold_CAD_Loft_Table.csv"):
    header = (
        f"Optimized Single-Inlet Regen Manifold Loft Table | "
        f"Optimal Inlet Wye Half-Angle: {opt_results['inlet_angle_deg']:.2f} deg | "
        f"Maldistribution: {opt_results['maldistribution_pct']:.3f}%\n"
        "Theta_deg,Centerline_X_mm,Centerline_Y_mm,Centerline_Z_mm,"
        "Inner_Radius_mm,Outer_Radius_mm,Wall_Thickness_mm,Flow_Area_mm2"
    )
    
    xyz_mm = cad_data["centerline_xyz"] * 1e3
    r_in_mm = cad_data["r_inner_full_m"] * 1e3
    r_out_mm = cad_data["r_outer_full_m"] * 1e3
    t_mm = cad_data["t_wall_full_m"] * 1e3
    area_mm2 = np.pi * r_in_mm**2
    
    table = np.column_stack((
        cad_data["theta_full_deg"],
        xyz_mm[:, 0],
        xyz_mm[:, 1],
        xyz_mm[:, 2],
        r_in_mm,
        r_out_mm,
        t_mm,
        area_mm2
    ))
    
    np.savetxt(filename, table, delimiter=",", header=header, comments="# ", fmt="%.5f")
    print(f"[Geometry] Saved 360-deg CAD Loft Table to: {filename}")