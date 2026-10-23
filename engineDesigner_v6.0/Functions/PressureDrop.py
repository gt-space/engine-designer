import numpy as np

#dP due to friction losses
def dPdarcyWeisbach(Re, rho, V, L, roughness, Dh):
    #darcy friction factor, Mileikovskyi Tkachenko
    A0 = -0.79638*np.log(roughness/Dh/8.208+7.3357/Re) # A0
    A1 = Re*roughness/Dh+9.3120665*A0 # A1
    f = ((8.128943 + A1)/(8.128943 * A0 - 0.86859209 * A1* np.log(A1 / (3.7099535 * Re))))**2
    
    #Apply Darcy Weisbach Equation
    dP = (L*f/Dh* rho * V**2 /2)
    return dP

def dPareaChange(Dh1, Dh2, rho, V):
    K = ((Dh1/Dh2)**2-1)**2
    dP = K*rho*V**2/2
    return dP

def dPfluidAcceleration(A_down, A_up, rho, mdot_ch):
    dP = mdot_ch**2 / (2*rho) * (
        1/A_down**2 - 1/A_up**2
    )
    return dP

