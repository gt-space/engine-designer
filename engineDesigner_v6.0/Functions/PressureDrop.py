import numpy as np

def darcyWeishbachFriction(Re, roughness, Dh):
#darcy friction factor, Mileikovskyi Tkachenko
    A0 = -0.79638*np.log(roughness/Dh/8.208+7.3357/Re) # A0
    A1 = Re*roughness/Dh+9.3120665*A0 # A1
    return ((8.128943 + A1)/(8.128943 * A0 - 0.86859209 * A1* np.log(A1 / (3.7099535 * Re))))**2

#dP due to friction losses
def dPdarcyWeisbach(Re, rho, V, L, roughness, Dh):
    #Apply Darcy Weisbach Equation
    f=darcyWeishbachFriction(Re,roughness,Dh)
    dP = (L*f/Dh* rho * V**2 /2)
    return dP

# Minor Losses (Bends, Tees, Contractions, Channel Entry/Exit via K-factor)
def dPminorLoss(K, rho, V):
    return K * 0.5 * rho * V**2

# Sudden Expansion / Contraction (Using True Area Ratio, referenced to V_upstream or V_downstream)
def dPareaChange(A_up, A_down, rho, V_ref, ref_is_downstream=True):
    if A_down > A_up:
        # Borda-Carnot Sudden Expansion (referenced to upstream higher velocity V_up)
        V_up = V_ref * (A_down / A_up) if ref_is_downstream else V_ref
        K_exp = (1.0 - A_up / A_down)**2
        return K_exp * 0.5 * rho * V_up**2
    elif A_down < A_up:
        # Sudden Contraction (Idelchik, referenced to downstream higher velocity V_down)
        V_down = V_ref if ref_is_downstream else V_ref * (A_up / A_down)
        K_con = 0.5 * (1.0 - A_down / A_up)**0.75
        return K_con * 0.5 * rho * V_down**2
    return 0.0

# Reversible Bernoulli Area + Thermal Density Acceleration dP
def dPfluidAcceleration(A_down, A_up, rho_down, rho_up, mdot_ch):
    rho_avg = 0.5 * (rho_up + rho_down)
    A_avg = 0.5 * (A_up + A_down)
    
    # Reversible Bernoulli static pressure change due to area change
    dP_area = (mdot_ch**2 / (2.0 * rho_avg)) * (1.0 / A_down**2 - 1.0 / A_up**2)
    # Momentum pressure drop due to thermal expansion (density decrease along channel)
    dP_thermal = (mdot_ch / A_avg)**2 * (1.0 / rho_down - 1.0 / rho_up)
    return dP_area + dP_thermal

# Distributing Manifold Friction & Static Pressure Regain 
def dPmanifoldSplit(Re_man, rho, V0, L_quad, roughness, Dh_man, beta=1.0):
    # 1/3 factor accounts for linear mass bleed-off along manifold channels
    dP_fric = (1.0 / 3.0) * dPdarcyWeisbach(Re_man, rho, V0, L_quad, roughness, Dh_man)
    dP_regain = beta * 0.5 * rho * V0**2
    return dP_fric - dP_regain



