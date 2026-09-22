
import math

def skin_friction(gas, x, local_radius):

    #local reynolds number
    x_safe = max(x, 1 * (10**-5))
    reynolds_local = (gas.rho * gas.core_velo * x_safe) / gas.mu

    #baseline skin friction
    #flat plate assumption
    skin_fric_o = 0.0592 * (reynolds_local ** -0.2)

    #dynamic turbulence
    d = 2 * local_radius
    reynolds_d = max((gas.rho * gas.core_velo * d) / gas.mu, 1)

    core_turb = 0.16 * (reynolds_d ** -0.125)
    core_shear = math.sqrt(skin_fric_o / 2)

    turb_intensity = math.sqrt((core_turb ** 2) + (core_shear ** 2))
    k_turb = 1 + (4 * turb_intensity)
    
    #baseline stanton num
    St_o = (skin_fric_o / 2) * (gas.Pr ** -0.667) * k_turb

    #mass transfer conductance
    g_m = gas.rho * gas.core_velo * St_o

    #shear stress from core gas
    shear = 0.5 * skin_fric_o * gas.rho * (gas.core_velo **2)

    return skin_fric_o, St_o, g_m, shear