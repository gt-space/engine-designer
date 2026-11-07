import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize

import Functions.Structural as st
import Solvers.ManifoldSolver as fn

#Interpolates a small set of radial control points (from 0 to pi radians) 
#onto all N_half cooling channel stations using a monotonic PCHIP spline.
def interpolateManifoldProfile(r_control_pts, N_half, shape="circular", aspect_ratio=1.0):

    n_ctrl = len(r_control_pts)
    theta_ctrl = np.linspace(0.0, np.pi, n_ctrl)
    theta_stations = np.linspace(0.0, np.pi, N_half)
    
    # PCHIP preserves monotonicity 
    spline = PchipInterpolator(theta_ctrl, r_control_pts)
    r_profile = np.maximum(spline(theta_stations), 1e-4)
    
    if shape == "rectangular":
        A_profile, Dh_profile = st.calcDuctSection("rectangular", 2.0 * r_profile, 2.0 * r_profile * aspect_ratio)
    else:
        A_profile, Dh_profile = st.calcDuctSection(shape, r_profile)
        
    return theta_stations, r_profile, A_profile, Dh_profile

#Makes Intial Guess to get constant velocity and match intial manifold velocity
def generateConstantVelocityInitialGuess(engine_params, n_ctrl=5, shape="circular", aspect_ratio=1.0):
    A_pipe = engine_params['A_feed_pipe']
    A_inlet_0 = 0.5 * A_pipe
    
    # Linear area taper from 100% at theta=0 down to 15% at theta=180 deg
    taper_fractions = np.linspace(1.0, 0.15, n_ctrl)
    A_ctrl = A_inlet_0 * taper_fractions
    
    if shape == "circular":
        r_ctrl_0 = np.sqrt(A_ctrl / np.pi)
    elif shape == "semicircular":
        r_ctrl_0 = np.sqrt(2.0 * A_ctrl / np.pi)
    elif shape == "rectangular":
        r_ctrl_0 = 0.5 * np.sqrt(A_ctrl / aspect_ratio)
    else:
        r_ctrl_0 = np.sqrt(A_ctrl / np.pi)
        
    inlet_angle_0 = 45.0  # Start at a 45-deg Y-inlet split
    return np.append(r_ctrl_0, inlet_angle_0)

#Main optimization routine
def optimizeManifold(engine_params, channel_geom, struct_params,
                     n_ctrl=5, shape="circular", aspect_ratio=1.0, max_iter=25):

    N_half_full = channel_geom["N_channels"] // 2

    # --- 1. Coarse Network Reduction Setup ---
    # Evaluating 45 FullFlow states for every gradient step is what caused the slow runtime.
    # Group into <= 15 representative segments during optimization, then re-evaluate full 45 at the end.
    N_half_opt = min(N_half_full, 15)
    scale_factor = N_half_full / N_half_opt

    fast_channel_geom = dict(channel_geom)
    fast_channel_geom["N_channels"] = N_half_opt * 2
    fast_channel_geom["A_ch_inlet"] = channel_geom["A_ch_inlet"] * scale_factor
    fast_channel_geom["A_channel_avg"] = channel_geom["A_channel_avg"] * scale_factor

    # --- 2. Build Initial Guess & Bounds ---
    x0 = generateConstantVelocityInitialGuess(engine_params, n_ctrl, shape, aspect_ratio)

    r_min = struct_params.get("r_min_mfg", 0.002)
    r_max = struct_params.get("max_radial_envelope", 0.04)
    angle_bounds = (20.0, 90.0)

    bounds = [(r_min, r_max) for _ in range(n_ctrl)] + [angle_bounds]

    history = {
        "iteration": [],
        "cost": [],
        "maldistribution_pct": [],
        "dP_total_bar": [],
        "dP_static_diff_bar": [],
        "inlet_angle_deg": [],
        "manifold_mass_kg": []
    }

    dP_scale = 2.0e5  # 2 bar reference pressure scale
    R_base_man = channel_geom["R_manifold"]

    def objectiveFunction(x):
        r_ctrl = x[:-1]
        inlet_angle = x[-1]

        # A. Interpolate manifold geometry profile on fast mesh
        theta_rad, r_prof, A_prof, Dh_prof = interpolateManifoldProfile(
            r_ctrl, N_half_opt, shape=shape, aspect_ratio=aspect_ratio
        )

        # B. Solve 1D FullFlow Fluid Network
        try:
            flow_res = fn.solveManifoldFlowNetwork(
                A_prof, Dh_prof, inlet_angle, engine_params, fast_channel_geom
            )
        except Exception:
            return 1e6

        # C. Structural Wall Thickness, Mass, and Outer Envelope Penalty
        P_peak = flow_res["P_inlet_required"]
        t_wall = st.calcStructuralWall(
            P_peak, r_prof,
            struct_params["yield_strength"],
            struct_params.get("safety_factor", 2.0),
            struct_params.get("t_min_mfg", 0.0015)
        )

        # Integrate full 360-deg ring (2x half-ring) with physical torus radius
        theta_half = np.linspace(0.0, np.pi, N_half_opt)
        _, metal_mass_half = st.calcRevolvedShellMass(
            theta_half, r_prof, t_wall, struct_params["rho_material"],
            path_type="circumferential", R_torus_base=R_base_man
        )
        metal_mass = 2.0 * metal_mass_half

        actual_protrusion = 2.0 * r_prof + t_wall
        env_penalty = st.calcEnvelopePenalty(
            actual_protrusion, struct_params["max_radial_envelope"], weight=500.0
        )

        # D. Cost Function
        maldist_pct = flow_res["maldistribution_pct"]
        cost_maldist = 50.0 * (maldist_pct / 1.0)**2

        dP_end_to_end = flow_res["dP_static_end_to_end"]
        cost_static_flatness = 20.0 * (dP_end_to_end / 1e4)**2

        dP_tot = flow_res["dP_total"]
        cost_dP = 5.0 * (dP_tot / dP_scale)**2

        growth_violations = np.maximum(0.0, np.diff(r_ctrl))
        cost_taper_smoothness = 100.0 * np.sum((growth_violations / r_min)**2)

        cost_mass = 1.0 * metal_mass

        total_cost = (
            cost_maldist
            + cost_static_flatness
            + cost_dP
            + cost_taper_smoothness
            + cost_mass
            + env_penalty
        )

        history["iteration"].append(len(history["iteration"]) + 1)
        history["cost"].append(total_cost)
        history["maldistribution_pct"].append(maldist_pct)
        history["dP_total_bar"].append(dP_tot / 1e5)
        history["dP_static_diff_bar"].append(dP_end_to_end / 1e5)
        history["inlet_angle_deg"].append(inlet_angle)
        history["manifold_mass_kg"].append(metal_mass)

        return total_cost

    # --- 3. Run Constrained Optimization ---
    res = minimize(
        objectiveFunction,
        x0,
        method="SLSQP",
        bounds=bounds,
        options={"maxiter": max_iter, "ftol": 1e-3, "eps": 1e-3, "disp": True}
    )

    # --- 4. Re-evaluate Final Optimal Geometry on the FULL Channel Network ---
    r_ctrl_opt = res.x[:-1]
    inlet_angle_opt = res.x[-1]

    theta_rad_full, r_opt_full, A_opt_full, Dh_opt_full = interpolateManifoldProfile(
        r_ctrl_opt, N_half_full, shape=shape, aspect_ratio=aspect_ratio
    )

    final_flow = fn.solveManifoldFlowNetwork(
        A_opt_full, Dh_opt_full, inlet_angle_opt, engine_params, channel_geom
    )

    t_wall_opt_full = st.calcStructuralWall(
        final_flow["P_inlet_required"], r_opt_full,
        struct_params["yield_strength"],
        struct_params.get("safety_factor", 2.0),
        struct_params.get("t_min_mfg", 0.0015)
    )

    return {
        "success": res.success,
        "message": res.message,
        "r_control_points": r_ctrl_opt,
        "inlet_angle_deg": inlet_angle_opt,
        "theta_stations_rad": theta_rad_full,
        "theta_stations_deg": np.degrees(theta_rad_full),
        "r_manifold_profile": r_opt_full,
        "A_manifold_profile": A_opt_full,
        "Dh_manifold_profile": Dh_opt_full,
        "t_wall_profile": t_wall_opt_full,
        "mdot_channels": final_flow["mdot_distribution"],
        "P_static_manifold": final_flow["P_static_manifold"],
        "maldistribution_pct": final_flow["maldistribution_pct"],
        "dP_total_Pa": final_flow["dP_total"],
        "P_inlet_required_Pa": final_flow["P_inlet_required"],
        "history": history
    }