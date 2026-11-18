import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize
from scipy.stats import norm, lognorm, qmc

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

    # Taper area down to 6% at 180 deg to match the linear bleed-off of mass flow
    taper_fractions = np.linspace(1.0, 0.06, n_ctrl)
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
    fast_geom["A_ch_inlet"] = channel_geom["A_ch_inlet"]       
    fast_geom["A_channel_avg"] = channel_geom["A_channel_avg"] 
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
                     max_maldistribution_pct=1.0, max_iter=50):

    N_half_full = channel_geom["N_channels"] // 2
    N_half_opt = min(N_half_full, 15)
    fast_channel_geom, f_mult = _prepare_fast_channel_geom(engine_params, channel_geom, N_half_opt)

    x0_si = generateConstantVelocityInitialGuess(engine_params, n_ctrl, shape, aspect_ratio)
    x0_scaled = np.append(x0_si[:-1] * 1e3, x0_si[-1] / 10.0)

    r_pipe_mm = np.sqrt(engine_params["A_feed_pipe"] / np.pi) * 1e3  
    r_min_mm = struct_params.get("r_min_mfg", 0.0018) * 1e3
    r_max_mm = min(struct_params.get("max_radial_envelope", 0.04) * 1e3, r_pipe_mm * 0.95)

    bounds = [(r_min_mm, r_max_mm) for _ in range(n_ctrl)] + [(3.0, 7.5)]

    def to_si(x_sc):
        return np.append(x_sc[:-1] * 1e-3, x_sc[-1] * 10.0)

    history = {
        "iteration": [], "cost": [], "maldistribution_pct": [],
        "dP_total_bar": [], "P_static_std_bar": [],
        "inlet_angle_deg": [], "manifold_mass_kg": []
    }

    def objective(x_sc):
        x_si = to_si(x_sc)
        try:
            flow_res, r_prof, _, _, t_wall, mass, _ = _evaluate_candidate(
                x_si, N_half_opt, engine_params, fast_channel_geom, struct_params, shape, aspect_ratio
            )
        except Exception:
            return 1e6

        maldist = flow_res["maldistribution_pct"]
        p_std = flow_res["P_static_std"]
        dP_tot = flow_res["dP_total"]

        # Smoothness regularizer on control point radii (prevents S-curve stair-stepping)
        r_mm = x_sc[:-1]
        d2_r = np.diff(r_mm, n=2)
        cost_smoothness = 0.15 * np.sum(d2_r ** 2)

        # Balance uniformity against total pressure drop and manifold mass
        cost_maldist = 20.0 * (maldist / max_maldistribution_pct) ** 2
        cost_p_std = 10.0 * (p_std / 5000.0) ** 2
        cost_dP = 8.0 * (dP_tot / 2.5e5)

        actual_protrusion = 2.0 * r_prof + t_wall
        env_pen = st.calcEnvelopePenalty(actual_protrusion, struct_params["max_radial_envelope"], weight=200.0)
        total_cost = cost_maldist + cost_p_std + cost_dP + cost_smoothness + 1.5 * mass + env_pen

        history["iteration"].append(len(history["iteration"]) + 1)
        history["cost"].append(total_cost)
        history["maldistribution_pct"].append(maldist)
        history["dP_total_bar"].append(dP_tot / 1e5)
        history["P_static_std_bar"].append(p_std / 1e5)
        history["inlet_angle_deg"].append(x_si[-1])
        history["manifold_mass_kg"].append(mass)

        return total_cost

    # Strictly decreasing taper (at least 0.3 mm drop per control segment)
    constraints = [
        {"type": "ineq", "fun": lambda x_sc: (x_sc[:-2] - x_sc[1:-1]) - 0.3}
    ]

    res_sc = minimize(
        objective, x0_scaled, method="SLSQP", bounds=bounds, constraints=constraints,
        options={"maxiter": max_iter, "ftol": 1e-4, "eps": 2e-2, "disp": True}
    )

    # Convert optimal vector back to SI meters & degrees before finalizing
    res_sc.x = to_si(res_sc.x)
    return _finalize_manifold_result(
        res_sc, engine_params, channel_geom, struct_params, f_mult, shape, aspect_ratio, history, "Deterministic"
    )

# 2. Optimization Under Uncertainity

# Defines Coefficient of Variation for Dh_channel, Roughness, rho, mu, plus RBDO targets
def getDefaultUncertaintyConfig(beta_target=2.0, pf_target=None):
    if pf_target is not None:
        # Convert probability of failure
        beta_val = float(-norm.ppf(pf_target))
    else:
        beta_val = float(beta_target)
        pf_target = float(norm.sf(beta_val))

    return {
        "cov_Dh_channel": 0.04,
        "cov_roughness": 0.20,
        "cov_rho": 0.02,
        "cov_mu": 0.08,
        "beta_target": beta_val,       # Reliability index target 
        "pf_target": pf_target,        # Allowable probability of failure P_f = Phi(-beta)
        "k_reliability": beta_val,     
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
# 2. Reliability-Based Design Optimization (RBDO) Under Uncertainty
def optimizeManifoldUnderUncertainty(engine_params, channel_geom, struct_params,
                                     uq_config=None, n_ctrl=5, shape="circular",
                                     aspect_ratio=1.0, max_maldistribution_pct=1.5,
                                     max_P_inlet_Pa=None, max_iter=30, x0_init=None):
    if uq_config is None:
        uq_config = getDefaultUncertaintyConfig()

    # Support either explicit beta_target or pf_target in uq_config
    if "pf_target" in uq_config and "beta_target" not in uq_config:
        beta_target = float(-norm.ppf(uq_config["pf_target"]))
    else:
        beta_target = float(uq_config.get("beta_target", uq_config.get("k_reliability", 2.0)))
    pf_target = float(norm.sf(beta_target))

    N_half_full = channel_geom["N_channels"] // 2
    N_half_opt = min(N_half_full, 15)
    fast_channel_geom, f_mult = _prepare_fast_channel_geom(engine_params, channel_geom, N_half_opt)

    sigma_pts, weights = _generate_unscented_sigma_points(engine_params, fast_channel_geom, uq_config)

    x0 = x0_init if x0_init is not None else generateConstantVelocityInitialGuess(
        engine_params, n_ctrl, shape, aspect_ratio
    )
    r_min = struct_params.get("r_min_mfg", 0.0025)
    r_max = struct_params.get("max_radial_envelope", 0.04)
    bounds = [(r_min, r_max) for _ in range(n_ctrl)] + [(20.0, 90.0)]

    dP_scale = 2.5e5
    P_std_scale = 1.5e4

    history = {
        "iteration": [], "cost": [],
        "maldistribution_mean_pct": [], "maldistribution_std_pct": [],
        "beta_maldist": [], "pf_maldist": [],
        "dP_mean_bar": [], "dP_std_bar": [],
        "inlet_angle_deg": [], "manifold_mass_kg": []
    }

    # Cache to avoid re-running the 9 UT sigma points when SLSQP queries objective & constraints at the same x
    _cache = {"x": None, "res": None}

    def _eval_ut_state(x):
        if _cache["x"] is not None and np.array_equal(x, _cache["x"]):
            return _cache["res"]

        r_ctrl = x[:-1]
        inlet_angle = x[-1]

        theta_rad, r_prof, A_prof, Dh_prof = interpolateManifoldProfile(
            r_ctrl, N_half_opt, shape=shape, aspect_ratio=aspect_ratio
        )

        dP_samples, p_std_samples, maldist_samples, P_in_samples = [], [], [], []

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
                P_in_samples.append(res_pt["P_inlet_required"])
            except Exception:
                state = {"valid": False}
                _cache["x"] = x.copy()
                _cache["res"] = state
                return state

        dP_arr = np.array(dP_samples)
        p_std_arr = np.array(p_std_samples)
        mal_arr = np.array(maldist_samples)
        P_in_arr = np.array(P_in_samples)

        # Unscented Transform Mean and Standard Deviation
        mean_dP = float(np.sum(weights * dP_arr))
        std_dP = float(np.sqrt(max(np.sum(weights * (dP_arr - mean_dP) ** 2), 1e-12)))

        mean_p_std = float(np.sum(weights * p_std_arr))
        std_p_std = float(np.sqrt(max(np.sum(weights * (p_std_arr - mean_p_std) ** 2), 1e-12)))

        mean_mal = float(np.sum(weights * mal_arr))
        std_mal = float(np.sqrt(max(np.sum(weights * (mal_arr - mean_mal) ** 2), 1e-8)))

        mean_P_in = float(np.sum(weights * P_in_arr))
        std_P_in = float(np.sqrt(max(np.sum(weights * (P_in_arr - mean_P_in) ** 2), 1e-12)))

        # Reliability Index (beta) and Probability of Failure (Pf) for maldistribution limit
        beta_mal = (max_maldistribution_pct - mean_mal) / std_mal
        pf_mal = float(norm.sf(beta_mal))

        # Structural wall & mass sized at beta_target upper-bound inlet pressure
        P_peak_rbdo = mean_P_in + beta_target * std_P_in
        t_wall = st.calcStructuralWall(
            P_peak_rbdo, r_prof,
            struct_params["yield_strength"],
            struct_params.get("safety_factor", 2.0),
            struct_params.get("t_min_mfg", 0.0015)
        )
        _, mass_half = st.calcRevolvedShellMass(
            theta_rad, r_prof, t_wall, struct_params["rho_material"],
            path_type="circumferential", R_torus_base=channel_geom["R_manifold"]
        )
        mass = 2.0 * mass_half

        state = {
            "valid": True,
            "r_prof": r_prof,
            "t_wall": t_wall,
            "mass": mass,
            "mean_dP": mean_dP, "std_dP": std_dP,
            "mean_p_std": mean_p_std, "std_p_std": std_p_std,
            "mean_mal": mean_mal, "std_mal": std_mal,
            "mean_P_in": mean_P_in, "std_P_in": std_P_in,
            "beta_mal": beta_mal, "pf_mal": pf_mal,
            "inlet_angle": inlet_angle
        }
        _cache["x"] = x.copy()
        _cache["res"] = state
        return state

    def rbdo_objective(x):
        state = _eval_ut_state(x)
        if not state["valid"]:
            return 1e6

        # RBDO Objective: Minimize expected (mean) & variance of performance metrics + structural mass
        # Feasibility/failure modes (maldistribution & envelope) are enforced via hard RBDO constraints
        robust_dP = state["mean_dP"] + 0.5 * state["std_dP"]
        robust_p_std = state["mean_p_std"] + 0.5 * state["std_p_std"]

        cost = (
            10.0 * (robust_dP / dP_scale)
            + 8.0 * (robust_p_std / P_std_scale) ** 2
            + 0.25 * state["mass"]
        )

        history["iteration"].append(len(history["iteration"]) + 1)
        history["cost"].append(cost)
        history["maldistribution_mean_pct"].append(state["mean_mal"])
        history["maldistribution_std_pct"].append(state["std_mal"])
        history["beta_maldist"].append(state["beta_mal"])
        history["pf_maldist"].append(state["pf_mal"])
        history["dP_mean_bar"].append(state["mean_dP"] / 1e5)
        history["dP_std_bar"].append(state["std_dP"] / 1e5)
        history["inlet_angle_deg"].append(state["inlet_angle"])
        history["manifold_mass_kg"].append(state["mass"])

        return cost

    # 1. Probabilistic Maldistribution Constraint:
    # Enforces mu_mal + beta_target * std_mal <= max_maldistribution_pct (Normalized to O(1) for SLSQP)
    def constraint_maldistribution_reliability(x):
        state = _eval_ut_state(x)
        if not state["valid"]:
            return -1.0
        g_mal = max_maldistribution_pct - (state["mean_mal"] + beta_target * state["std_mal"])
        return g_mal / max_maldistribution_pct

    # 2. Radial Envelope Constraint (with RBDO-sized wall thickness):
    # Enforces 2*r_prof + t_wall <= max_radial_envelope at all stations
    def constraint_radial_envelope(x):
        state = _eval_ut_state(x)
        if not state["valid"]:
            return -1.0
        max_env = struct_params["max_radial_envelope"]
        actual_protrusion = 2.0 * state["r_prof"] + state["t_wall"]
        return (max_env - np.max(actual_protrusion)) / max_env

    constraints = [
        # Monotonic taper constraint: r_1 >= r_2 >= ... >= r_n
        {"type": "ineq", "fun": lambda x: (x[:-2] - x[1:-1]) * 1e3},
        # Hard RBDO Maldistribution constraint (beta >= beta_target / Pf <= pf_target)
        {"type": "ineq", "fun": constraint_maldistribution_reliability},
        # Hard Radial Envelope constraint
        {"type": "ineq", "fun": constraint_radial_envelope},
    ]

    # 3. Optional Probabilistic Peak Inlet Pressure Constraint:
    # Enforces mu_P_in + beta_target * std_P_in <= max_P_inlet_Pa
    if max_P_inlet_Pa is not None:
        def constraint_inlet_pressure_reliability(x):
            state = _eval_ut_state(x)
            if not state["valid"]:
                return -1.0
            g_pin = max_P_inlet_Pa - (state["mean_P_in"] + beta_target * state["std_P_in"])
            return g_pin / dP_scale

        constraints.append({"type": "ineq", "fun": constraint_inlet_pressure_reliability})

    res = minimize(
        rbdo_objective, x0, method="SLSQP", bounds=bounds, constraints=constraints,
        options={"maxiter": max_iter, "ftol": 1e-4, "eps": 5e-4, "disp": True}
    )

    # Evaluate final UT reliability metrics at the converged optimum x_opt
    final_ut = _eval_ut_state(res.x)

    out = _finalize_manifold_result(
        res, engine_params, channel_geom, struct_params, f_mult, shape, aspect_ratio, history, "RBDO"
    )
    out["uq_config"] = uq_config
    out["beta_target"] = beta_target
    out["pf_target"] = pf_target
    if final_ut["valid"]:
        out["beta_maldist_UT"] = final_ut["beta_mal"]
        out["pf_maldist_UT"] = final_ut["pf_mal"]
        out["maldistribution_mean_UT_pct"] = final_ut["mean_mal"]
        out["maldistribution_std_UT_pct"] = final_ut["std_mal"]
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
                                uq_config=None, max_maldistribution_pct=1.5,
                                max_P_inlet_Pa=None, n_samples=100, seed=42):

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

    P_in_arr = np.array(P_in_list)
    dP_arr = np.array(dP_list)
    mal_arr = np.array(mal_list)
    p_std_arr = np.array(p_std_list)

    # 1. Empirical failure probability (fraction of LHS samples exceeding limit)
    pf_mal_empirical = float(np.mean(mal_arr > max_maldistribution_pct))
    # Clip between [1e-6, 1 - 1e-6] so norm.ppf remains finite when 0 samples fail
    beta_mal_empirical = float(-norm.ppf(np.clip(pf_mal_empirical, 1e-6, 1.0 - 1e-6)))

    # 2. First-order moment-based reliability index from LHS sample mean & std
    mal_mean = float(np.mean(mal_arr))
    mal_std = float(max(np.std(mal_arr, ddof=1), 1e-8))
    beta_mal_moment = (max_maldistribution_pct - mal_mean) / mal_std
    pf_mal_moment = float(norm.sf(beta_mal_moment))

    mc_results = {
        "P_inlet_Pa": P_in_arr,
        "dP_total_Pa": dP_arr,
        "maldistribution_pct": mal_arr,
        "P_static_std_Pa": p_std_arr,
        "maldistribution_mean_pct": mal_mean,
        "maldistribution_std_pct": mal_std,
        "pf_maldist_empirical": pf_mal_empirical,
        "beta_maldist_empirical": beta_mal_empirical,
        "pf_maldist_moment": pf_mal_moment,
        "beta_maldist_moment": beta_mal_moment,
    }

    if max_P_inlet_Pa is not None:
        pf_pin_emp = float(np.mean(P_in_arr > max_P_inlet_Pa))
        pin_mean = float(np.mean(P_in_arr))
        pin_std = float(max(np.std(P_in_arr, ddof=1), 1e-8))
        beta_pin_mom = (max_P_inlet_Pa - pin_mean) / pin_std
        mc_results["pf_P_inlet_empirical"] = pf_pin_emp
        mc_results["beta_P_inlet_moment"] = beta_pin_mom
        mc_results["pf_P_inlet_moment"] = float(norm.sf(beta_pin_mom))

    return mc_results