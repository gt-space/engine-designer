import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize

import Functions.Structural as st
import Functions.PressureDrop as pd
import Solvers.ManifoldSolver as fn

#Interpolates a small set of radial control points (from 0 to pi radians) 
#onto all N_half cooling channel stations using a monotonic PCHIP spline.
def interpolateManifoldProfile(r_control_pts, N_half, shape="circular", aspect_ratio=1.0):
    n_ctrl = len(r_control_pts)
    theta_ctrl = np.linspace(0.0, np.pi, n_ctrl)
    theta_stations = np.linspace(0.0, np.pi, N_half)

    spline = PchipInterpolator(theta_ctrl, r_control_pts)
    r_profile = np.maximum(spline(theta_stations), 1e-4)

    if shape == "rectangular":
        A_profile, Dh_profile = st.calcDuctSection("rectangular", 2.0 * r_profile, 2.0 * r_profile * aspect_ratio)
    else:
        A_profile, Dh_profile = st.calcDuctSection(shape, r_profile)

    return theta_stations, r_profile, A_profile, Dh_profile

#Makes Intial Guess to get constant velocity and match intial manifold velocity
def generateConstantVelocityInitialGuess(engine_params, n_ctrl=5, shape="circular", aspect_ratio=1.0):
    A_pipe = engine_params["A_feed_pipe"]
    A_inlet_0 = 0.55 * A_pipe

    # Taper area from 100% at inlet (0 deg) to 20% at 180 deg
    taper_fractions = np.linspace(1.0, 0.20, n_ctrl)
    A_ctrl = A_inlet_0 * taper_fractions

    if shape == "circular":
        r_ctrl_0 = np.sqrt(A_ctrl / np.pi)
    elif shape == "semicircular":
        r_ctrl_0 = np.sqrt(2.0 * A_ctrl / np.pi)
    elif shape == "rectangular":
        r_ctrl_0 = 0.5 * np.sqrt(A_ctrl / aspect_ratio)
    else:
        r_ctrl_0 = np.sqrt(A_ctrl / np.pi)

    inlet_angle_0 = 45.0
    return np.append(r_ctrl_0, inlet_angle_0)

def _prepare_fast_channel_geom(engine_params, channel_geom, N_half_opt):
    N_half_full = channel_geom["N_channels"] // 2
    scale_factor = N_half_full / N_half_opt

    # Calibrate f_Save_multiplier so DarcyWeisbach matches the coupled ChannelSolver nominal dP
    Dh_nom = channel_geom["Dh_channel"]
    V_nom = (0.5 * engine_params["mdot_total"] / N_half_full) / (
        engine_params["rho"] * channel_geom["A_channel_avg"]
    )
    Re_nom = engine_params["rho"] * V_nom * Dh_nom / engine_params["mu"]
    f_base_nom = float(pd.darcyWeishbachFriction(Re_nom, engine_params["roughness"], Dh_nom))
    f_mult = channel_geom.get("f_channel", f_base_nom) / max(f_base_nom, 1e-5)

    fast_geom = dict(channel_geom)
    fast_geom["N_channels"] = N_half_opt * 2
    fast_geom["scale_factor"] = scale_factor
    fast_geom["Dh_channel_nominal"] = Dh_nom
    fast_geom["f_Save_multiplier"] = f_mult
    fast_geom["A_ch_inlet"] = channel_geom["A_ch_inlet"] * scale_factor
    fast_geom["A_channel_avg"] = channel_geom["A_channel_avg"] * scale_factor
    return fast_geom, f_mult

def _evaluate_candidate(x, N_half_eval, engine_params, channel_geom, struct_params,
                        shape="circular", aspect_ratio=1.0):
    r_ctrl = x[:-1]
    inlet_angle = x[-1]

    theta_rad, r_prof, A_prof, Dh_prof = interpolateManifoldProfile(
        r_ctrl, N_half_eval, shape=shape, aspect_ratio=aspect_ratio
    )
    flow_res = fn.solveManifoldFlowNetwork(
        A_prof, Dh_prof, inlet_angle, engine_params, channel_geom
    )

    P_peak = flow_res["P_inlet_required"]
    t_wall = st.calcStructuralWall(
        P_peak, r_prof,
        struct_params["yield_strength"],
        struct_params.get("safety_factor", 2.0),
        struct_params.get("t_min_mfg", 0.0015)
    )

    _, metal_mass_half = st.calcRevolvedShellMass(
        theta_rad, r_prof, t_wall, struct_params["rho_material"],
        path_type="circumferential", R_torus_base=channel_geom["R_manifold"]
    )
    metal_mass = 2.0 * metal_mass_half

    # Reference dynamic pressure at manifold inlet
    V_in = (0.5 * engine_params["mdot_total"]) / (engine_params["rho"] * A_prof[0])
    q_dyn_in = max(0.5 * engine_params["rho"] * V_in**2, 1e3)

    return flow_res, r_prof, A_prof, Dh_prof, t_wall, metal_mass, q_dyn_in

# 1.Deteministic Manifold Optimization

#Deterministic optimization using nominal fluid properties and dimensions:
#Minimizes overall pressure requirement (dP_total) while maintaining a relatively uniform static pressure distribution
def optimizeManifold(engine_params, channel_geom, struct_params,
                     n_ctrl=5, shape="circular", aspect_ratio=1.0,
                     max_maldistribution_pct=1.5, max_iter=35):

    N_half_full = channel_geom["N_channels"] // 2
    N_half_opt = min(N_half_full, 15)
    fast_channel_geom, f_mult = _prepare_fast_channel_geom(engine_params, channel_geom, N_half_opt)

    x0 = generateConstantVelocityInitialGuess(engine_params, n_ctrl, shape, aspect_ratio)
    r_min = struct_params.get("r_min_mfg", 0.0025)
    r_max = struct_params.get("max_radial_envelope", 0.04)
    bounds = [(r_min, r_max) for _ in range(n_ctrl)] + [(20.0, 90.0)]

    dP_scale = 2.5e5  # 2.5 bar reference scale
    P_std_scale = 1.5e4  # 0.15 bar static pressure spatial uniformity scale

    history = {
        "iteration": [], "cost": [], "maldistribution_pct": [],
        "dP_total_bar": [], "P_static_std_bar": [],
        "inlet_angle_deg": [], "manifold_mass_kg": []
    }

    def evaluate_metrics(x):
        flow_res, r_prof, _, _, t_wall, mass, _ = _evaluate_candidate(
            x, N_half_opt, engine_params, fast_channel_geom, struct_params, shape, aspect_ratio
        )
        return flow_res, r_prof, t_wall, mass

    def objective(x):
        try:
            flow_res, r_prof, t_wall, mass = evaluate_metrics(x)
        except Exception:
            return 1e6

        # 1. Minimize overall inlet pressure requirement (dP_total = P_inlet - P_exit)
        dP_tot = flow_res["dP_total"]
        cost_pressure_req = 10.0 * (dP_tot / dP_scale)

        # 2. Maintain relatively uniform static pressure distribution via spatial std dev
        # (Does NOT force P_0 == P_end; allows deceleration pressure recovery)
        p_std = flow_res["P_static_std"]
        cost_static_uniformity = 8.0 * (p_std / P_std_scale) ** 2

        # 3. Channel flow maldistribution penalty + soft barrier above threshold
        maldist = flow_res["maldistribution_pct"]
        cost_maldist = 15.0 * (maldist / max_maldistribution_pct) ** 2

        # 4. Mass and envelope
        actual_protrusion = 2.0 * r_prof + t_wall
        env_pen = st.calcEnvelopePenalty(actual_protrusion, struct_params["max_radial_envelope"], weight=200.0)
        total_cost = cost_pressure_req + cost_static_uniformity + cost_maldist + 0.25 * mass + env_pen

        history["iteration"].append(len(history["iteration"]) + 1)
        history["cost"].append(total_cost)
        history["maldistribution_pct"].append(maldist)
        history["dP_total_bar"].append(dP_tot / 1e5)
        history["P_static_std_bar"].append(p_std / 1e5)
        history["inlet_angle_deg"].append(x[-1])
        history["manifold_mass_kg"].append(mass)

        return total_cost

    # Monotonic taper constraint: r_1 >= r_2 >= r_3 >= r_4 >= r_5
    constraints = [
        {"type": "ineq", "fun": lambda x: (x[:-2] - x[1:-1]) * 1e3}
    ]

    res = minimize(
        objective, x0, method="SLSQP", bounds=bounds, constraints=constraints,
        options={"maxiter": max_iter, "ftol": 1e-4, "eps": 5e-4, "disp": True}
    )

    return _finalize_manifold_result(res, engine_params, channel_geom, struct_params, f_mult, shape, aspect_ratio, history, "Deterministic")

# 2. Optimization Under Uncertainity

#Defines Coefficient of Variation for Dh_channel,Roughness, rho, mu
def getDefaultUncertaintyConfig():

    return {
        "cov_Dh_channel": 0.04,
        "cov_roughness": 0.20,
        "cov_rho": 0.02,
        "cov_mu": 0.08,
        "k_reliability": 2.0,  # Mean + 2*Std (~97.7% confidence bound)
    }

#Generates 2*n + 1 = 9 Unscented Transform (UT) sigma points in the 4D uncertain space
def _generate_unscented_sigma_points(engine_params, channel_geom, uq_config):

    means = np.array([
        channel_geom["Dh_channel"],
        engine_params["roughness"],
        engine_params["rho"],
        engine_params["mu"]
    ], dtype=float)

    stds = means * np.array([
        uq_config["cov_Dh_channel"],
        uq_config["cov_roughness"],
        uq_config["cov_rho"],
        uq_config["cov_mu"]
    ], dtype=float)

    n_dim = 4
    kappa = 2.0  # Standard Gaussian kurtosis scaling (n + kappa = 6)
    scale = np.sqrt(n_dim + kappa)

    sigma_pts = [means.copy()]
    weights = [kappa / (n_dim + kappa)]

    for i in range(n_dim):
        e_i = np.zeros(n_dim)
        e_i[i] = stds[i] * scale
        sigma_pts.append(np.maximum(means + e_i, 1e-7))
        sigma_pts.append(np.maximum(means - e_i, 1e-7))
        weights.extend([0.5 / (n_dim + kappa), 0.5 / (n_dim + kappa)])

    return np.array(sigma_pts), np.array(weights)

#Optimization Under Uncertainty
def optimizeManifoldUnderUncertainty(engine_params, channel_geom, struct_params,
                                     uq_config=None, n_ctrl=5, shape="circular",
                                     aspect_ratio=1.0, max_maldistribution_pct=1.5,
                                     max_iter=30, x0_init=None):
    if uq_config is None:
        uq_config = getDefaultUncertaintyConfig()

    N_half_full = channel_geom["N_channels"] // 2
    N_half_opt = min(N_half_full, 15)
    fast_channel_geom, f_mult = _prepare_fast_channel_geom(engine_params, channel_geom, N_half_opt)

    sigma_pts, weights = _generate_unscented_sigma_points(engine_params, fast_channel_geom, uq_config)
    k_rel = uq_config.get("k_reliability", 2.0)

    x0 = x0_init if x0_init is not None else generateConstantVelocityInitialGuess(engine_params, n_ctrl, shape, aspect_ratio)
    r_min = struct_params.get("r_min_mfg", 0.0025)
    r_max = struct_params.get("max_radial_envelope", 0.04)
    bounds = [(r_min, r_max) for _ in range(n_ctrl)] + [(20.0, 90.0)]

    dP_scale = 2.5e5
    P_std_scale = 1.5e4

    history = {
        "iteration": [], "cost": [], "maldistribution_mean_pct": [],
        "maldistribution_std_pct": [], "dP_mean_bar": [], "dP_std_bar": [],
        "inlet_angle_deg": []
    }

    def robust_objective(x):
        r_ctrl = x[:-1]
        inlet_angle = x[-1]

        theta_rad, r_prof, A_prof, Dh_prof = interpolateManifoldProfile(
            r_ctrl, N_half_opt, shape=shape, aspect_ratio=aspect_ratio
        )

        dP_samples = []
        p_std_samples = []
        maldist_samples = []

        for pt in sigma_pts:
            Dh_val, eps_val, rho_val, mu_val = pt

            ep_pert = dict(engine_params)
            ep_pert["roughness"] = float(eps_val)
            ep_pert["rho"] = float(rho_val)
            ep_pert["mu"] = float(mu_val)

            cg_pert = dict(fast_channel_geom)
            cg_pert["Dh_channel"] = float(Dh_val)

            try:
                res_pt = fn.solveManifoldFlowNetwork(A_prof, Dh_prof, inlet_angle, ep_pert, cg_pert)
                dP_samples.append(res_pt["dP_total"])
                p_std_samples.append(res_pt["P_static_std"])
                maldist_samples.append(res_pt["maldistribution_pct"])
            except Exception:
                return 1e6

        dP_arr = np.array(dP_samples)
        p_std_arr = np.array(p_std_samples)
        mal_arr = np.array(maldist_samples)

        # Compute Unscented Mean and Standard Deviation
        mean_dP = float(np.sum(weights * dP_arr))
        std_dP = float(np.sqrt(max(np.sum(weights * (dP_arr - mean_dP) ** 2), 0.0)))

        mean_p_std = float(np.sum(weights * p_std_arr))
        std_p_std = float(np.sqrt(max(np.sum(weights * (p_std_arr - mean_p_std) ** 2), 0.0)))

        mean_mal = float(np.sum(weights * mal_arr))
        std_mal = float(np.sqrt(max(np.sum(weights * (mal_arr - mean_mal) ** 2), 0.0)))

        # Structural mass at 97.7% upper-bound inlet pressure
        P_peak_robust = engine_params["P_collector_exit"] + mean_dP + k_rel * std_dP
        t_wall = st.calcStructuralWall(
            P_peak_robust, r_prof,
            struct_params["yield_strength"],
            struct_params.get("safety_factor", 2.0),
            struct_params.get("t_min_mfg", 0.0015)
        )
        _, mass_half = st.calcRevolvedShellMass(
            theta_rad, r_prof, t_wall, struct_params["rho_material"],
            path_type="circumferential", R_torus_base=channel_geom["R_manifold"]
        )
        mass = 2.0 * mass_half

        # Robust Cost: Mean + k * Std
        robust_dP = mean_dP + k_rel * std_dP
        robust_p_std = mean_p_std + k_rel * std_p_std
        robust_mal = mean_mal + k_rel * std_mal

        cost = (
            10.0 * (robust_dP / dP_scale)
            + 8.0 * (robust_p_std / P_std_scale) ** 2
            + 18.0 * (robust_mal / max_maldistribution_pct) ** 2
            + 0.25 * mass
        )

        history["iteration"].append(len(history["iteration"]) + 1)
        history["cost"].append(cost)
        history["maldistribution_mean_pct"].append(mean_mal)
        history["maldistribution_std_pct"].append(std_mal)
        history["dP_mean_bar"].append(mean_dP / 1e5)
        history["dP_std_bar"].append(std_dP / 1e5)
        history["inlet_angle_deg"].append(inlet_angle)

        return cost

    constraints = [
        {"type": "ineq", "fun": lambda x: (x[:-2] - x[1:-1]) * 1e3}
    ]

    res = minimize(
        robust_objective, x0, method="SLSQP", bounds=bounds, constraints=constraints,
        options={"maxiter": max_iter, "ftol": 1e-4, "eps": 5e-4, "disp": True}
    )

    out = _finalize_manifold_result(res, engine_params, channel_geom, struct_params, f_mult, shape, aspect_ratio, history, "OUU")
    out["uq_config"] = uq_config
    return out


def _finalize_manifold_result(res, engine_params, channel_geom, struct_params, f_mult, shape, aspect_ratio, history, mode_label):
    r_ctrl_opt = res.x[:-1]
    inlet_angle_opt = res.x[-1]
    N_half_full = channel_geom["N_channels"] // 2

    full_channel_geom = dict(channel_geom)
    full_channel_geom["scale_factor"] = 1.0
    full_channel_geom["Dh_channel_nominal"] = channel_geom["Dh_channel"]
    full_channel_geom["f_Save_multiplier"] = f_mult

    theta_rad_full, r_opt_full, A_opt_full, Dh_opt_full = interpolateManifoldProfile(
        r_ctrl_opt, N_half_full, shape=shape, aspect_ratio=aspect_ratio
    )

    final_flow = fn.solveManifoldFlowNetwork(
        A_opt_full, Dh_opt_full, inlet_angle_opt, engine_params, full_channel_geom
    )

    t_wall_opt_full = st.calcStructuralWall(
        final_flow["P_inlet_required"], r_opt_full,
        struct_params["yield_strength"],
        struct_params.get("safety_factor", 2.0),
        struct_params.get("t_min_mfg", 0.0015)
    )

    return {
        "mode": mode_label,
        "success": res.success,
        "message": res.message,
        "x_opt": res.x,
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
        "P_static_std_Pa": final_flow["P_static_std"],
        "maldistribution_pct": final_flow["maldistribution_pct"],
        "dP_total_Pa": final_flow["dP_total"],
        "P_inlet_required_Pa": final_flow["P_inlet_required"],
        "history": history
    }



# 3. Latin Hypercube Monte Carlo Validation

#Evaluates a fixed manifold geometry across N Latin Hypercube samples of (Dh_channel, roughness, rho, mu) to quantify output distributions.
def evaluateDesignMonteCarloLHS(opt_result, engine_params, channel_geom,
                                uq_config=None, n_samples=100, seed=42):

    if uq_config is None:
        uq_config = getDefaultUncertaintyConfig()

    N_half_full = channel_geom["N_channels"] // 2
    N_half_mc = min(N_half_full, 15)
    fast_geom, _ = _prepare_fast_channel_geom(engine_params, channel_geom, N_half_mc)

    r_ctrl = opt_result["r_control_points"]
    inlet_angle = opt_result["inlet_angle_deg"]
    _, _, A_prof, Dh_prof = interpolateManifoldProfile(r_ctrl, N_half_mc)

    sampler = qmc.LatinHypercube(d=4, seed=seed)
    u_samples = sampler.random(n=n_samples)

    from scipy.stats import norm, lognorm
    Dh_mean = channel_geom["Dh_channel"]
    Dh_samples = norm.ppf(u_samples[:, 0], loc=Dh_mean, scale=Dh_mean * uq_config["cov_Dh_channel"])

    # Lognormal for roughness so roughness is strictly positive
    eps_mean = engine_params["roughness"]
    s_ln = np.sqrt(np.log(1.0 + uq_config["cov_roughness"] ** 2))
    eps_samples = lognorm.ppf(u_samples[:, 1], s=s_ln, scale=eps_mean * np.exp(-0.5 * s_ln ** 2))

    rho_mean = engine_params["rho"]
    rho_samples = norm.ppf(u_samples[:, 2], loc=rho_mean, scale=rho_mean * uq_config["cov_rho"])

    mu_mean = engine_params["mu"]
    mu_samples = norm.ppf(u_samples[:, 3], loc=mu_mean, scale=mu_mean * uq_config["cov_mu"])

    P_in_list, dP_list, mal_list, p_std_list = [], [], [], []

    for i in range(n_samples):
        ep_i = dict(engine_params)
        ep_i["roughness"] = float(max(eps_samples[i], 1e-7))
        ep_i["rho"] = float(max(rho_samples[i], 100.0))
        ep_i["mu"] = float(max(mu_samples[i], 1e-5))

        cg_i = dict(fast_geom)
        cg_i["Dh_channel"] = float(max(Dh_samples[i], 1e-4))

        res_i = fn.solveManifoldFlowNetwork(A_prof, Dh_prof, inlet_angle, ep_i, cg_i)
        P_in_list.append(res_i["P_inlet_required"])
        dP_list.append(res_i["dP_total"])
        mal_list.append(res_i["maldistribution_pct"])
        p_std_list.append(res_i["P_static_std"])

    return {
        "P_inlet_Pa": np.array(P_in_list),
        "dP_total_Pa": np.array(dP_list),
        "maldistribution_pct": np.array(mal_list),
        "P_static_std_Pa": np.array(p_std_list),
    }