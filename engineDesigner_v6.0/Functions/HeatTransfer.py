
from EngineParameters import Pc
from Functions.Correlations import calculate_nusselt
import numpy as np
from Solvers.HotGasFlowSolver import getGasPropertiesCEA
from Geometry.Contour import getCurrentEngineParameters

Pc_Pa = Pc * 6894.75729


_CACHED_GAS_PROPS = None
_CACHED_GEOM_PARAMS = None


def _get_cached_bartz_constants():
    global _CACHED_GAS_PROPS, _CACHED_GEOM_PARAMS
    if _CACHED_GAS_PROPS is None:
        _CACHED_GAS_PROPS = getGasPropertiesCEA()
    if _CACHED_GEOM_PARAMS is None:
        _CACHED_GEOM_PARAMS = getCurrentEngineParameters()
    return _CACHED_GAS_PROPS, _CACHED_GEOM_PARAMS


def clear_heat_transfer_cache():
    """Call this if contour or chamber pressure changes mid-script."""
    global _CACHED_GAS_PROPS, _CACHED_GEOM_PARAMS
    _CACHED_GAS_PROPS = None
    _CACHED_GEOM_PARAMS = None


# Coolant-side convection coefficient
def coolant_htc(Re, Pr, k_cool, roughness, Dh, method="gnielinski", **kwargs):
    if method == "gnielinski":
        Nu = calculate_nusselt("gnielinski", Re=Re, Pr=Pr, roughness=roughness, Dh=Dh)
    else:
        Nu = calculate_nusselt(method, Re=Re, Pr=Pr, **kwargs)

    hc = Nu * k_cool / Dh
    return hc, Nu


# Gas-side Bartz convection coefficient
def bartz_hg(A_local, M_local, Twg):
    gas, geometry_params = _get_cached_bartz_constants()

    T0 = gas["T0"]
    gamma = gas["gamma"]
    cp_g = gas["cp_g"]
    mu_g = gas["mu_g"]
    Pr_g = gas["Pr"]
    Cstar_actual = gas["Cstar_actual"]

    R_t = geometry_params["R_t"] * 0.0254   # in -> m
    Rtc_Rt = geometry_params["Rtc_Rt"]
    Rtc = Rtc_Rt * R_t
    At = np.pi * R_t**2

    # Find Bartz Correction Factor
    mach_term = 1.0 + 0.5 * (gamma - 1.0) * M_local**2
    sigma = 1.0 / ((0.5 * (Twg / T0) * mach_term + 0.5)**0.68 * mach_term**0.12)

    # Find Bartz Gas-side Convection Coefficient
    hg = (
        0.026 / (2.0 * R_t)**0.2
        * (mu_g**0.2 * cp_g / Pr_g**0.6)
        * (Pc_Pa / Cstar_actual)**0.8
        * ((2.0 * R_t) / Rtc)**0.1
        * (At / A_local)**0.9
        * sigma
    )
    return hg


def ievlev_hg(B, tau_lambda, D_bar, P_cc, d_th, S_K_T, Pr):
    term1 = B * (tau_lambda / (D_bar**1.82))
    term2 = (P_cc**0.85) / (d_th**0.15)
    term3 = S_K_T / (Pr**0.58)
    return term1 * term2 * term3


def fin_efficiency(hc, kw, h_ch, w_fin):
    m = np.sqrt(2.0 * hc / (kw * w_fin))
    Lc = h_ch + 0.5 * w_fin
    N_fin = np.tanh(m * Lc) / (m * Lc)
    return N_fin