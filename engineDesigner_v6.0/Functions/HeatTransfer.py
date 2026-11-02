from EngineParameters import Pc
from Functions.Correlations import calculate_nusselt
import numpy as np
from Solvers.HotGasFlowSolver import getGasPropertiesCEA

# Convert chamber pressure to Pascals once
Pc_Pa = Pc * 6894.75729

#Coolant-side convection coefficient
def coolant_htc(Re, Pr, k_cool, roughness, Dh, method="gnielinski", **kwargs):
    if method == "gnielinski":
        Nu = calculate_nusselt("gnielinski", Re=Re, Pr=Pr, roughness=roughness, Dh=Dh)
    else:
        Nu = calculate_nusselt(method, Re=Re, Pr=Pr, **kwargs)

    hc = Nu * k_cool / Dh
    return hc, Nu

#Gas-side Bartz convection coefficient 
def bartz_hg(A_local, M_local,Twg):
    
    #Get gas properties
    gas = getGasPropertiesCEA()
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


    #Find Bartz Correction Factor
    mach_term = (1+(gamma-1)/2*M_local**2)
    sigma = 1/((0.5* (Twg / T0)* mach_term+ 0.5)**(0.68) * mach_term**(0.12))
    #Find Bartz Gas-side Convection Coefficient
    hg = (
        0.026 / (2*R_t)**0.2
        * (mu_g**0.2 * cp_g / Pr_g**0.6)
        * (Pc_Pa / Cstar_actual)**0.8
        * ((2*R_t) / (Rtc))**0.1
        * (At / A_local)**0.9
        * sigma
    )
    return hg

def ievlev_hg(B, tau_lambda, D_bar, P_cc, d_th, S_K_T, Pr):
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

def fin_efficiency(hc, kw, h_ch, w_fin):
    m = np.sqrt(2*hc/(kw*w_fin))
    Lc = h_ch+w_fin/2
    N_fin = np.tanh(m*Lc)/(m*Lc)
    return N_fin