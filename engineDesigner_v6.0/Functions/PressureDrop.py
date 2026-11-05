import numpy as np

def darcyWeishbachFriction(Re, roughness, Dh):
    Re = np.maximum(Re, 1e-5)
    if np.isscalar(Re) and Re < 2300.0:
        return 64.0 / Re
    A0 = -0.79638 * np.log(roughness / Dh / 8.208 + 7.3357 / Re)
    A1 = Re * roughness / Dh + 9.3120665 * A0
    f_turb = ((8.128943 + A1) / (8.128943 * A0 - 0.86859209 * A1 * np.log(A1 / (3.7099535 * Re))))**2
    if not np.isscalar(Re):
        return np.where(Re < 2300.0, 64.0 / Re, f_turb)
    return f_turb

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

#dP at the single inlet split for manifold
def dPinletSplit(mdot_total, rho, A_pipe, A_man_inlet, inlet_angle_deg=90.0):

    V_pipe = mdot_total / (rho * A_pipe)
    V_man_0 = (0.5 * mdot_total) / (rho * A_man_inlet)
    
    alpha_rad = np.radians(np.clip(inlet_angle_deg, 15.0, 90.0))
    vel_ratio = V_man_0 / np.maximum(V_pipe, 1e-6)
    
    # Idelchik dividing flow loss coefficient K referenced to upstream pipe dynamic head
    K_split = 1.0 - (2.0 * vel_ratio * np.cos(alpha_rad)) + vel_ratio**2
    
    # Add sudden area change loss if manifold total inlet area (2 * A_man_inlet) != A_pipe
    dP_turn = dPminorLoss(K_split, rho, V_pipe)
    dP_area = dPareaChange(A_pipe, 2.0 * A_man_inlet, rho, V_man_0, ref_is_downstream=True)
    
    dP_total_split = dP_turn + dP_area
    
    # Static pressure drop includes reversible Bernoulli velocity head change
    dP_bernoulli = 0.5 * rho * (V_man_0**2 - V_pipe**2)
    dP_static_split = dP_total_split + dP_bernoulli
    
    return dP_total_split, dP_static_split

#dP along one discrete angular segment of the manifold  
def dPmanifoldSegment(mdot_up, mdot_down, A_up, A_down, Dh_up, Dh_down, dL, rho, mu, roughness, gamma=0.75):
    V_up = mdot_up / (rho * np.maximum(A_up, 1e-9))
    V_down = mdot_down / (rho * np.maximum(A_down, 1e-9))
    
    V_avg = 0.5 * (V_up + V_down)
    Dh_avg = 0.5 * (Dh_up + Dh_down)
    
    Re_seg = (rho * V_avg * Dh_avg) / mu
    
    # 1. Irreversible wall friction drop along segment dL
    dP_fric_loss = dPdarcyWeisbach(Re_seg, rho, V_avg, dL, roughness, Dh_avg)
    
    # 2. Reversible momentum change (Deceleration -> Static Pressure Regain)
    # If V_down < V_up, dP_momentum_regain is POSITIVE (pressure rises)
    # If V_down > V_up (over-tapered), use standard Bernoulli (gamma = 1.0)
    recovery_factor = gamma if V_down <= V_up else 1.0
    dP_momentum_regain = recovery_factor * 0.5 * rho * (V_up**2 - V_down**2)
    
    # Net change in static pressure from station i to i+1
    dP_static_change = dP_momentum_regain - dP_fric_loss
    
    return dP_static_change, dP_fric_loss

#dP for 90 turn from manifold to cooling channel
def dPchannelBranch(mdot_man_local, mdot_ch, A_man_local, A_ch_inlet, rho, K_entry_base=0.5):

    V_man = mdot_man_local / (rho * np.maximum(A_man_local, 1e-9))
    V_ch = mdot_ch / (rho * np.maximum(A_ch_inlet, 1e-9))
    
    # 1. Turning loss due to cross-flow momentum in manifold 
    vel_ratio_sq = (V_man / np.maximum(V_ch, 1e-6))**2
    K_turn = 1.0 + vel_ratio_sq
    dP_turn = dPminorLoss(K_turn, rho, V_ch)
    
    # 2. Area contraction loss from manifold cross-section to channel slot
    dP_contract = dPareaChange(A_man_local, A_ch_inlet, rho, V_ch, ref_is_downstream=True)
    
    # 3. Reversible Bernoulli static-to-dynamic conversion entering the channel
    dP_bernoulli_entry = 0.5 * rho * V_ch**2 * (1.0 + K_entry_base)
    
    return dP_turn + dP_contract + dP_bernoulli_entry

