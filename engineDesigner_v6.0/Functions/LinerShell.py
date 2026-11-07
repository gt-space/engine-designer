import numpy as np
from EngineParameters import load_matprotlib_material
from Functions.Structural import queryMatPropArray, calcStructuralWall



# 1. Inner Liner Thickness Equations (Structural vs. Thermal)

#Calculates required minimum inner liner thickness based on structural mechanics
def calcLinerStructuralThickness(P_cool, P_gas, r_inner, w_ch, material, Tw_avg,
                                 safety_factor=2.0, t_min_mfg=0.75e-3):
    mat = load_matprotlib_material(material)
    sigma_y = queryMatPropArray(mat, "yield_strength", Tw_avg, fallback_val=250e6)
    sigma_allow = np.maximum(sigma_y / safety_factor, 1e6)

    dP_slot = np.abs(P_cool - P_gas)
    t_bending = w_ch * np.sqrt((0.5 * dP_slot) / sigma_allow)
    t_shear = (dP_slot * w_ch) / (2.0 * 0.577 * sigma_allow)

    dP_burst = np.maximum(P_gas - P_cool, 0.0)
    t_hoop = (dP_burst * r_inner) / sigma_allow

    return np.maximum.reduce([t_bending, t_shear, t_hoop, np.full_like(w_ch, t_min_mfg)])

#Calculates required inner liner thickness based on coupled thermo-structural properties
def calcLinerThermalThickness(P_cool, P_gas, w_ch, q_flux, Twg, Twc, material,
                              safety_factor=1.5, t_min_mfg=0.75e-3):

    mat = load_matprotlib_material(material)
    Tw_avg = 0.5 * (Twg + Twc)

    sigma_y_hot = queryMatPropArray(mat, "yield_strength", Twg, fallback_val=150e6)
    E_hot = queryMatPropArray(mat, "elastic_modulus", Twg, fallback_val=100e9)
    kw_avg = queryMatPropArray(mat, "thermal_conductivity", Tw_avg, fallback_val=300.0)
    nu = queryMatPropArray(mat, "poisson_ratio", Tw_avg, fallback_val=0.33)
    cte = queryMatPropArray(mat, "cte", Tw_avg, fallback_val=16.5e-6)

    T_melt = queryMatPropArray(mat, "melting_point", Tw_avg, fallback_val=1350.0)
    try:
        T_max_meta = mat.get_meta("max_operating_temp")
        T_limit = float(T_max_meta) if T_max_meta is not None else 0.75 * T_melt
    except Exception:
        T_limit = 0.75 * T_melt

    dT_wall = np.maximum(Twg - Twc, 1.0)
    thermal_strain_factor = (E_hot * cte * dT_wall) / (
        2.0 * np.maximum(1.0 - nu, 0.1) * np.maximum(sigma_y_hot, 1e6)
    )
    dP_slot = np.abs(P_cool - P_gas)

    sigma_allow_hot = np.maximum(sigma_y_hot / safety_factor, 1e6)
    t_thermal_req = (
        w_ch
        * np.sqrt((0.75 * dP_slot) / sigma_allow_hot)
        * np.sqrt(1.0 + np.clip(thermal_strain_factor, 0.0, 2.0))
    )
    t_thermal_req = np.maximum(t_thermal_req, t_min_mfg)

    # 1D Fourier Law conduction ceiling before hot-gas wall overheats beyond T_limit
    q_safe = np.maximum(q_flux, 1e3)
    dT_headroom = np.maximum(T_limit - Twc, 10.0)
    t_max_fourier = (kw_avg * dT_headroom) / q_safe

    return t_thermal_req, t_max_fourier



# 2. Outer Closeout Shell Equations (Structural vs. Thermal)

#Calculates required outer closeout shell thickness after cooling channels based on structural loads
def calcShellStructuralThickness(P_cool, r_shell_inner, w_ch, shell_material, T_shell,
                                 safety_factor=2.0, t_min_mfg=1.5e-3):

    mat = load_matprotlib_material(shell_material)
    sigma_y = queryMatPropArray(mat, "yield_strength", T_shell, fallback_val=800e6)
    sigma_allow = np.maximum(sigma_y / safety_factor, 1e6)

    t_vessel = calcStructuralWall(
        P_cool, r_shell_inner, sigma_y,
        safety_factor=safety_factor, t_min_mfg=t_min_mfg
    )
    t_slot_bridge = w_ch * np.sqrt((0.5 * P_cool) / sigma_allow)

    return np.maximum.reduce([t_vessel, t_slot_bridge, np.full_like(r_shell_inner, t_min_mfg)])

#Calculates required outer closeout shell thickness based on material props
def calcShellThermalThickness(P_cool, r_shell_inner, t_liner, h_ch, f_fin,
                              Twg, Twc, T_cool, liner_material, shell_material,
                              safety_factor=1.5, t_min_mfg=1.5e-3, T_ref=293.15):
    
    mat_liner = load_matprotlib_material(liner_material)
    mat_shell = load_matprotlib_material(shell_material)

    T_liner_avg = 0.5 * (Twg + Twc)
    T_shell = np.maximum(T_cool, T_ref)

    E_liner = queryMatPropArray(mat_liner, "elastic_modulus", T_liner_avg, fallback_val=110e9)
    cte_liner = queryMatPropArray(mat_liner, "cte", T_liner_avg, fallback_val=17.0e-6)

    sigma_y_shell = queryMatPropArray(mat_shell, "yield_strength", T_shell, fallback_val=900e6)
    E_shell = queryMatPropArray(mat_shell, "elastic_modulus", T_shell, fallback_val=200e9)
    cte_shell = queryMatPropArray(mat_shell, "cte", T_shell, fallback_val=13.0e-6)

    eps_liner = cte_liner * np.maximum(T_liner_avg - T_ref, 0.0)
    eps_shell = cte_shell * np.maximum(T_shell - T_ref, 0.0)
    d_eps_th = np.abs(eps_liner - eps_shell)

    t_liner_eff = t_liner + f_fin * h_ch
    sigma_allow_shell = np.maximum(sigma_y_shell / safety_factor, 1e6)

    stiffness_ratio = E_liner / np.maximum(E_shell, 1e9)
    # Coaxial compliance factor (prevents over-predicting thickness when E_liner == E_shell)
    thermal_line_load = 0.25 * (E_liner * t_liner_eff * d_eps_th) / (1.0 + stiffness_ratio)
    pressure_axial_line_load = 0.5 * P_cool * r_shell_inner

    t_shell_thermal = (thermal_line_load + pressure_axial_line_load) / sigma_allow_shell
    return np.maximum(t_shell_thermal, t_min_mfg)