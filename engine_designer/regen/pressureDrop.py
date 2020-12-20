# Scirpt to calculate the pressure lost at a station along the engine
from .colebrook import colebrook

def pressureDrop(Re, D_hyd, vel, rho, length):
    k = 0.005 * (10 ** (-3)) #Surface roughness, guess based on our manufacuting capability, mm => m

    # Solve for Darcy Friction Factor Approximation using Reynolds number
    if Re <= 2320:
        fric = 64 / Re # Laminar flow approximation
    else:
        fric = colebrook(Re, D_hyd, k) # Colebrook-White Equation solution

    dP = fric * (length/D_hyd) * rho * ((vel**2)/2) #Solve for final dP across passage segment (Pa)
    return dP
