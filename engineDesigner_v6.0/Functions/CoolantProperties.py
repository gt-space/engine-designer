import numpy as np
from CoolProp.CoolProp import PropsSI

# Options: "RP1_NIST" (True MIL-DTL-25576 RP-1 correlations) or "n-Dodecane" (CoolProp surrogate)
COOLANT_MODEL = "RP1_NIST"
fluid = "n-Dodecane"

# Reference conditions for NIST RP-1 enthalpy integration
_T_REF = 298.15       # K
_P_REF = 101325.0     # Pa


def _rp1_density(T, P):
    dT = T - _T_REF
    # 1 atm density polynomial (270 K to 650 K)
    rho_0 = 806.2 - 0.732 * dT - 2.45e-4 * (dT**2)
    rho_0 = np.maximum(rho_0, 350.0)

    # Isothermal bulk modulus B(T) [Pa] for high-pressure compression (P up to 35 MPa)
    B_T = np.maximum(1.32e9 * np.exp(-0.0042 * dT), 1.5e8)
    dP = np.maximum(P - _P_REF, 0.0)
    return rho_0 / (1.0 - 0.0894 * np.log(1.0 + dP / B_T))

#Isobaric specific heat capacity [J/(kg*K)] of liquid/supercritical RP-1.
def _rp1_cp(T):

    # Linear/quadratic fit from 270 K (1985 J/kg-K) to 650 K (3310 J/kg-K)
    return 1985.0 + 3.65 * (T - _T_REF) + 4.5e-4 * ((T - _T_REF)**2)

#  Dynamic viscosity [Pa*s] of RP-1 with Vogel-Fulcher-Tammann (VFT) temperature dependence and Barus pressure correction.
def _rp1_viscosity(T, P):
  
    T_safe = np.maximum(T, 240.0)
    # VFT equation calibrated to NIST RP-1 (1.68 mPa*s at 298.15 K)
    mu_0 = 3.12e-5 * np.exp(612.0 / (T_safe - 144.5))
    # High-pressure viscosity enhancement factor
    pressure_factor = np.exp(1.1e-8 * np.maximum(P - _P_REF, 0.0))
    return mu_0 * pressure_factor

#Thermal conductivity [W/(m*K)] of RP-1 with pressure enhancement.
def _rp1_conductivity(T, P):
    dT = T - _T_REF
    k_0 = np.maximum(0.1145 - 1.48e-4 * dT + 2.5e-8 * (dT**2), 0.055)
    # Pressure correction (~4% increase per 10 MPa)
    k_P = k_0 * (1.0 + 3.8e-9 * np.maximum(P - _P_REF, 0.0))
    return k_P


def get_coolant_h(T, P):
    if COOLANT_MODEL == "n-Dodecane":
        return PropsSI("Hmass", "T", T, "P", P, fluid)

    # Analytical integral of cp(T) dT + incompressible Poynting pressure work (P - P_ref) / rho
    dT = T - _T_REF
    h_thermal = 1985.0 * dT + 0.5 * 3.65 * (dT**2) + (4.5e-4 / 3.0) * (dT**3)
    rho_avg = _rp1_density(0.5 * (T + _T_REF), P)
    h_pressure = (P - _P_REF) / rho_avg
    return h_thermal + h_pressure


def get_coolant_props(P, h):
    if COOLANT_MODEL == "n-Dodecane":
        T = PropsSI("T", "P", P, "Hmass", h, fluid)
        rho = PropsSI("Dmass", "P", P, "Hmass", h, fluid)
        mu = PropsSI("V", "P", P, "Hmass", h, fluid)
        k = PropsSI("L", "P", P, "Hmass", h, fluid)
        cp = PropsSI("Cpmass", "P", P, "Hmass", h, fluid)
        Pr = cp * mu / k
        return {"T": T, "rho": rho, "mu": mu, "k": k, "cp": cp, "Pr": Pr}

    # Fast 3-step Newton-Raphson inversion of h(T, P) -> T
    T_est = _T_REF + (h - (P - _P_REF) / 806.2) / 2100.0
    for _ in range(4):
        f_val = get_coolant_h(T_est, P) - h
        df_dT = _rp1_cp(T_est)
        T_est -= f_val / df_dT

    T = float(np.clip(T_est, 240.0, 900.0))
    rho = float(_rp1_density(T, P))
    cp = float(_rp1_cp(T))
    mu = float(_rp1_viscosity(T, P))
    k = float(_rp1_conductivity(T, P))
    Pr = cp * mu / k

    return {
        "T": T,
        "rho": rho,
        "mu": mu,
        "k": k,
        "cp": cp,
        "Pr": Pr,
    }
