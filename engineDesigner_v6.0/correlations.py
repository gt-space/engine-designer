import math
import numpy as np

# ==========================================
# IMPORTING CRITICAL DATABASES & LIBRARIES
# ==========================================
import matprotlib as mp           # Material Database
from rocketcea.cea_obj import CEA_Obj  # Rocket CEA Database

try:
    import CoolProp.CoolProp as CP 
except ImportError:
    pass

# ==========================================
# 1. HOT-GAS SIDE HEAT TRANSFER (COMBUSTION)
# ==========================================

def heat_flux_temperature(h_g, T_gaw, T_gw):
    """Hot-Gas Side Convective Heat Flux using Temperature."""
    return h_g * (T_gaw - T_gw)

def heat_flux_enthalpy(h_g, c_p_gn, i_gaw, i_gw):
    """Hot-Gas Side Convective Heat Flux using Enthalpy."""
    return (h_g / c_p_gn) * (i_gaw - i_gw)

def adiabatic_wall_enthalpy(i_gs, Pr_gx, i_g0):
    """Adiabatic wall enthalpy (Bartz and Eckert)."""
    return i_gs + (Pr_gx**(1/3)) * (i_g0 - i_gs)

def bartz_sigma_correction(T_hw, T_s, gamma, M):
    """Calculates the sigma correction factor for the Bartz correlation."""
    term1 = 0.5 * (T_hw / T_s) * (1 + ((gamma - 1) / 2) * M**2) + 0.5
    term2 = 1 + ((gamma - 1) / 2) * M**2
    return 1.0 / ((term1**0.68) * (term2**0.12))

def bartz_correlation(D_star, mu_0, cp_0, Pr_0, p_c, g, c_star, r_c, A_star, A, sigma):
    """
    Full Bartz Correlation for Combustion Gases.
    Notes:
    - Highly sensitive to hot gas transport properties.
    - Best for rapid, preliminary design and initial throat sizing.
    - Most accurate at and near the nozzle throat.
    - Known under-predictor in the combustion chamber (often 40-50% off actuals).
    - Ignores radiation and combustion instabilities.
    """
    term1 = 0.026 / (D_star**0.2)
    term2 = ((mu_0**0.2) * cp_0) / (Pr_0**0.6)
    term3 = ((p_c * g) / c_star)**0.8
    term4 = (D_star / r_c)**0.1
    term5 = (A_star / A)**0.9
    return term1 * term2 * term3 * term4 * term5 * sigma

def ievlev_correlation(B, tau_lambda, D_bar, P_cc, d_th, S_K_T, Pr):
    """
    Ievlev semi-analytical convective heat flux.
    Notes:
    - A known over-predictor (originally designed for pipes, treats the engine like a pipe).
    - Very similar to Bartz (< 10% difference in the throat area).
    - Highly accurate at the throat but decreases significantly away from it.
    """
    # Note: tau_lambda is divided by D_bar**1.82 to match the math logic properly
    term1 = B * (tau_lambda / (D_bar**1.82))
    term2 = (P_cc**0.85) / (d_th**0.15)
    term3 = S_K_T / (Pr**0.58)
    return term1 * term2 * term3

# ==========================================
# 2. ADVANCED GAS PROPERTIES (KINETIC THEORY)
# ==========================================
# Introduced for Coupled Heat Transfer Calculations 

def eucken_thermal_conductivity(mu_i, W_i, c_pi):
    """
    Calculates the thermal conductivity coefficient of a gas component using Eucken's formula[cite: 1126, 1601].
    Useful for varying properties in fluid-structure coupled domains.
    """
    R_0 = 8.314462618 # Universal gas constant in J/(mol K)
    return (mu_i * R_0 / W_i) * (c_pi * W_i / R_0 + 1.25)

def lennard_jones_viscosity(W_i, T, sigma_i, epsilon_i):
    """
    Calculates the dynamic viscosity coefficient using the Lennard-Jones formula[cite: 1124, 1605].
    sigma_i = Lennard-Jones collision diameter (Angstroms)[cite: 1607].
    epsilon_i = Lennard-Jones well depth[cite: 1614].
    """
    k_B = 1.380649e-23 # Boltzmann constant [cite: 1612, 1614]
    T_star = (k_B * T) / epsilon_i # Reduced temperature [cite: 1611, 1613]
    
    # Collision integral calculation [cite: 1608, 1610]
    omega_mu = 1.147 * (T_star**-0.145) + (T_star + 0.5)**-2 
    
    # Dynamic Viscosity output [cite: 1605, 1607]
    return 2.6693e-6 * math.sqrt(W_i * T) / ((sigma_i**2) * omega_mu)


# ==========================================
# 3. COOLANT LIQUIDS (LAMINAR FLOW - FORCED & MIXED)
# ==========================================

def hausen_laminar(D, L, Re, Pr, mu, mu_w):
    """
    Hausen correlation for laminar forced convection.
    Notes: Underpredicts compared to experimental data. Fails to accurately predict transition.
    """
    numerator = 0.0668 * (D / L) * Re * Pr
    denominator = 1 + 0.04 * ((D / L) * Re * Pr)**(2/3)
    return 3.66 + (numerator / denominator) * ((mu / mu_w)**0.14)

def kays_laminar(D, L, Re, Pr, mu, mu_w):
    """
    Kays correlation for laminar forced convection.
    Notes: Valid only for small temperature differences between wall and bulk.
    """
    term1 = (0.36 * Re * Pr) / (L / D)
    term2 = math.log((L / D) / (0.0011 * Re * Pr) + 1)
    return 4.36 + term1 * term2 * ((mu / mu_w)**0.14)

def colburn_mixed(Re, Pr, D, L, mu, mu_w, Gr):
    """
    Colburn correlation for laminar mixed convection.
    Notes: Most consistent for Laminar region Re < 2300 (deviations 6.5% to 19.8%). 
    Overpredicts significantly in the entrance region (initial 100 diameters).
    """
    term1 = 1.5 * ((Re * Pr * D / L)**(1/3))
    buoyancy_term = 1 + 0.015 * (Gr**(1/3))
    return term1 * ((mu / mu_w)**(1/3)) * buoyancy_term

def mori_mixed(Gr, Pr, Re):
    """
    Mori et al. correlation for laminar mixed convection.
    Notes: Had high errors (10.2% to 27.8% deviation) in coiled/spiral pipe tests.
    """
    base_term = (Gr * Pr * Re)**0.2
    return 0.6 * base_term * (1 + (1.8 / base_term))

def eubank_proctor_mixed(Pr, Re, D, L, Gr, mu, mu_w):
    """
    Eubank and Proctor correlation for laminar mixed convection.
    Notes: Improvement to Sieder and Tate incorporating free convection. 
    Aims for ~25% deviation, but experimental checks hit up to 43% off.
    """
    term1 = Pr * Re * D / L
    term2 = 0.04 * (Gr * Pr * D / L)**0.75
    return 1.75 * ((term1 + term2)**(1/3)) * ((mu / mu_w)**0.14)


# ==========================================
# 4. COOLANT LIQUIDS (TURBULENT FLOW)
# ==========================================

def dittus_boelter_turbulent(Re, Pr):
    """
    Standard Dittus-Boelter correlation.
    Notes: Foundational correlation (originally for automobile coolant).
    Overpredicts heat transfer data (deviations 2.0% to 21.7%).
    """
    return 0.023 * (Re**0.8) * (Pr**0.4)

def sieder_tate_turbulent(Re, Pr, mu, mu_w):
    """
    Sieder-Tate correlation for turbulent forced convection.
    Notes: Works best for turbulent/upper-transition regions. 
    Underpredicts at very high Reynolds numbers. Tightest deviation range: 2.5% to 13.9%.
    """
    return 0.023 * (Re**0.8) * (Pr**(1/3)) * ((mu / mu_w)**0.14)

def sleicher_rouse_turbulent(Re, Pr):
    """Sleicher and Rouse correlation. Over-predictor (deviations 3.1% to 21.9%)."""
    a = 0.88 - (0.24 / (4 + Pr))
    b = 0.333 + 0.5 * math.exp(-0.6 * Pr)
    return 5 + 0.015 * (Re**a) * (Pr**b)

def petukhov_popov_turbulent(Re, Pr, mu, mu_w, heating=True):
    """
    Petukhov and Popov correlation for variable physical properties.
    Notes: Set heating=True (n=0.11) when T_wall > T_bulk; False (n=0.25) for cooling.
    """
    f = (1.82 * math.log10(Re) - 1.64)**-2
    a = 1.07 + 12.7 * (Pr**(2/3) - 1) * (f / 8)**0.5
    n = 0.11 if heating else 0.25 
    return ((Re * Pr * (f / 8)) / a) * ((mu / mu_w)**n)

def mcadams_turbulent_a(Re, Pr, L, D):
    """McAdams variant A. Overpredictor (deviations 3.0% to 21.9%)."""
    return 0.023 * (Re**0.8) * (Pr**0.4) * (1 + 1.4 / (L / D))

def mcadams_turbulent_b(Re, Pr, L, D):
    """McAdams variant B. Overpredictor (deviations 1.9% to 21.7%)."""
    return 0.023 * (Re**0.8) * (Pr**0.4) * (1 + ((L / D)**-0.7))

def gnielinski_turbulent(Re, Pr):
    """
    Gnielinski correlation.
    Notes: Recommended for Re > 2300. Overpredictor, but consistent within 10%.
    """
    f = (1.82 * math.log10(Re) - 1.64)**-2
    numerator = (f / 8) * (Re - 1000) * Pr
    denominator = 1 + 12.7 * ((f / 8)**0.5) * (Pr**(2/3) - 1)
    return numerator / denominator

def friend_metzner_turbulent(Re, Pr, mu, mu_w):
    """
    Friend and Metzner correlation.
    Notes: Consistently overpredicted data. Poor results in lower transition region.
    """
    return 0.015 * (Re**0.83) * (Pr**0.42) * ((mu / mu_w)**0.14)

def churchill_turbulent_a(Re, Pr):
    """
    Churchill correlation variant A.
    Notes: Underpredicts data heavily (40-80%) for Re < 4000. Depends on inlet geometry.
    """
    f_05 = 1 / (2.21 * math.log(Re / 7))
    numerator = 0.079 * f_05 * Re * Pr
    denominator = (1 + Pr**0.8)**(5/6)
    return 6.3 + (numerator / denominator)


# ==========================================
# 5. SPECIFIC PROPELLANTS & GASES
# ==========================================

def rp1_shell_correlation(Re_cs, Pr_cs):
    """Shell correlation specifically tailored for RP-1 liquid coolant."""
    return 0.255 * (Re_cs**0.582) * (Pr_cs**0.554)

def rp1_rocketdyne_correlation(Re_cx, Pr_cx):
    """Rocketdyne correlation specifically tailored for RP-1 liquid coolant."""
    return 0.0055 * (Re_cx**0.95) * (Pr_cx**0.4)

def liquid_hydrogen_htc(C_c, Re, Pr, gamma, T_w, T_s):
    """Cooling channel heat transfer coefficient for Liquid Hydrogen."""
    psi = 1 + gamma * (T_w - T_s)
    Nu_r = psi**-0.55
    return Nu_r * C_c * (Re**0.8) * (Pr**0.4)

def taylor_gaseous(Re, Pr, T_w, T_b, D_h, x):
    """Taylor correlation for gaseous coolants (H2, N2, Air, CO2)."""
    exponent = -(0.57 - 1.59 * (D_h / x))
    return 0.023 * (Re**0.8) * (Pr**0.4) * ((T_w / T_b)**exponent)

def mccarthy_wolf_hydrogen(Re, Pr, T_w, T_b):
    """McCarthy and Wolf correlation specifically for High Pressure Hydrogen."""
    return 0.025 * (Re**0.8) * (Pr**0.4) * ((T_w / T_b)**-0.55)



# ==========================================
# 6. THE CORRELATION DISPATCHER
# ==========================================

CORRELATIONS = {
    # Laminar Forced & Mixed
    "hausen": hausen_laminar,
    "kays": kays_laminar,
    "colburn": colburn_mixed,
    "mori": mori_mixed,
    "eubank_proctor": eubank_proctor_mixed,
    
    # Turbulent Forced
    "sieder_tate": sieder_tate_turbulent,
    "dittus_boelter": dittus_boelter_turbulent,
    "sleicher_rouse": sleicher_rouse_turbulent,
    "petukhov_popov": petukhov_popov_turbulent,
    "mcadams_a": mcadams_turbulent_a,
    "mcadams_b": mcadams_turbulent_b,
    "gnielinski": gnielinski_turbulent,
    "friend_metzner": friend_metzner_turbulent,
    "churchill_a": churchill_turbulent_a,
    
    # Specific Propellants
    "rp1_shell": rp1_shell_correlation,
    "rp1_rocketdyne": rp1_rocketdyne_correlation,
    "liquid_hydrogen": liquid_hydrogen_htc,
    
    # Gaseous
    "taylor": taylor_gaseous,
    "mccarthy_wolf": mccarthy_wolf_hydrogen
}

def calculate_nusselt(method_name, **kwargs):
    """
    Master dispatcher function to easily swap out correlations without rewriting your loops.
    
    Example Usage:
    Nu = calculate_nusselt("dittus_boelter", Re=50000, Pr=0.7)
    Nu = calculate_nusselt("sieder_tate", Re=50000, Pr=0.7, mu=1.2e-3, mu_w=1.0e-3)
    Nu = calculate_nusselt("colburn", Re=1500, Pr=0.7, D=0.01, L=1.0, mu=1.2e-3, mu_w=1.0e-3, Gr=1e5)
    """
    method_name = method_name.lower()
    if method_name not in CORRELATIONS:
        raise ValueError(f"Correlation '{method_name}' not found. Available options: {list(CORRELATIONS.keys())}")
    
    # Calls the mapped function with the provided arguments
    return CORRELATIONS[method_name](**kwargs)