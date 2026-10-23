from Geometry.Contour import getCurrentEngineParameters
from Solvers.HotGasFlowSolver import getGasPropertiesCEA
from EngineParameters import Pc
import numpy as np
from CoolProp.CoolProp import PropsSI

#Get gas properties
gas = getGasPropertiesCEA()
T0 = gas["T0"]
gamma = gas["gamma"]
cp_g = gas["cp_g"]
mu_g = gas["mu_g"]
Pr_g = gas["Pr"]
Cstar_actual = gas["Cstar_actual"]

#Get engine geometry
geometry_params = getCurrentEngineParameters()
R_t = geometry_params["R_t"] * 0.0254   # in -> m
Rtc_Rt = geometry_params["Rtc_Rt"]
Rtc = Rtc_Rt * R_t
At = np.pi * R_t**2

#Convert chamber pressure to pascals
Pc_Pa = Pc*6894.75729


#Coolant-side convection coefficient
def coolant_htc(Re, Pr, k_cool, roughness, Dh):

    #don't accept channel designs that have a reynolds number < 3000
    #caucses gnielinksi correlation to fail
    if Re <= 3000:
        raise ValueError(f"Reynolds number too low for Gnielinski correlation: Re = {Re:.1f}")
    
    #darcy friction factor, Mileikovskyi Tkachenko
    A0 = -0.79638*np.log(roughness/Dh/8.208+7.3357/Re) # A0
    A1 = Re*roughness/Dh+9.3120665*A0 # A1
    f = ((8.128943 + A1)/(8.128943 * A0 - 0.86859209 * A1* np.log(A1 / (3.7099535 * Re))))**2
    
    #Nusselts Number
    Nu = ((f/8)*(Re-1000)*Pr/(1+ 12.7 * np.sqrt(f / 8)*(Pr**(2/3)-1))) #Gnielinski

    #Convection Coefficient
    hc = Nu * k_cool / Dh
    return hc, Nu

#gas side Bartz convection coefficient 
def bartz_hg(A_local, M_local,Twg):

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

def fin_efficiency(hc, kw, h_ch, w_fin):
    m = np.sqrt(2*hc/(kw*w_fin))
    Lc = h_ch+w_fin/2
    N_fin = np.tanh(m*Lc)/(m*Lc)
    return N_fin