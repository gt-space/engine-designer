# Scirpt to calculate the pressure lost at a station along the engine
from colebrook import colebrook


def pressureDrop(CSA, perim, mDot_chan, rho_c, len, viscK_c):
    # pressureDrop - Calculate pressure drop in coolant passage

    epsilon = 0.005 * (10 ** (-3)) #Surface roughness, guess based on our manufacuting capability, mm
    d_hyd = 4 * CSA / perim #Equivalent hydraulic diameter of coolant passage
    K = epsilon / d_hyd
    w_c = mDot_chan/(rho_c * CSA) #Coolant Velocity

    Re_c = (w_c * d_hyd)/viscK_c #Reynolds number for coolant in passage (characteristic length chosen to be hydraulic diameter

    # Solve for Darcy Friction Factor Approximation using Reynolds number
    if Re_c <= 2320:
        fric = 64 / Re_c #Laminar flow approximation
    else:
        fric = colebrook(Re_c, K) #Colebrook-White Equation solution (via online function)

    dP = fric * (len/d_hyd) * rho_c * ((w_c**2)/2) #Solve for final dP across passage segment (Pa)

    return dP
