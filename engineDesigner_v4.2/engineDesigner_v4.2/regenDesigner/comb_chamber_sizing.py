from contourDesigner.design import Engine
import numpy as np
from rocketcea.cea_obj import CEA_Obj


def comb_chamber_sizing(engine, axial_vel):
    r_0 = 0.00005 #Initial droplet radius
    c_p = engine.engineProps[0, 14]*1000
    mu = engine.engineProps[0, 17]/1000 * 0.1 #mPoise -> Poise -> kg/ms
    cond = 0.1 * engine.engineProps[0, 18]
    gam = engine.engineProps[0, 15]

    rho_f= 800
    rho_o = 1141
    rho = engine.engineProps[0, 10]

    #Fix this
    h_fg_f = 251 #kJ/kg
    c_p_f = 2010 #J/kgK
    fuel_temp = 300 #K
    B_f = engine.engineProps[0, 14] * (engine.engineProps[0, 9] - fuel_temp) / h_fg_f

    h_fg_o = 6.82 / (32/1000) #kJ/mol * (mol/kg)
    c_p_o = 54
    lox_temp = 100
    B_o = engine.engineProps[0, 14] * (engine.engineProps[0, 9] - lox_temp) / h_fg_o

    Gamma = ((gam+1)/2)**((gam+1)/(2*(gam-1)))
    a = np.sqrt(gam * 8314.5/engine.engineProps[0, 13] * engine.engineProps[0, 9])
    phi_f = 9/2 * engine.engineProps[0, 19] / B_f
    phi_o = 9/2 * engine.engineProps[0, 19] / B_o

    B_f = 10
    phi_f = 0.5

    length_f = r_0**2 * (axial_vel/a + 3/Gamma * 1/con_rat * phi_f/10) * c_p * rho_f/cond * a/np.log(1+B_f) * 1/(2 + phi_f)
    #length_f = r_0 ** 2 * (axial_vel / a + 3 / Gamma * 1 / con_rat * phi_f / 10) * a / (2*10**-7) * 1 / (2 + phi_f)
    print(length_f)






#Injector parameters
#Fuel Side
Cd_f = 0.85
A_f = 3.249*10**-5
dP_f = 360000
rho_f = 800
mdot_f = Cd_f * A_f * np.sqrt(2*dP_f*rho_f)
vel_f = Cd_f * np.sqrt(2*dP_f / rho_f)
#LOX side
Cd_o = 0.7638
A_o = 5.7448139*10**-5
dP_o = 324000
rho_o = 1141
mdot_o = Cd_o * A_o * np.sqrt(2*dP_o*rho_o)
vel_o = Cd_o * np.sqrt(2*dP_o / rho_o)

axial_vel = vel_f * mdot_f/(mdot_o + mdot_f)
## GENERAL INPUT PARAMETERS ##
# For basic engine configuration and adjustment
thrust = 3500 #781 * 4.44822 # Thrust [lbf to N]
P_c = 18 #261 * 0.0689476 # Chamber pressure at injector face [psi to bar]
P_e = 8 / 14.5038 # Desired exit presure for expansion ratio [bar] (if in doubt put ambient)
con_rat = 6 # Contraction ratio
L_star = 1.05 # Characteristic length [m] (for recommended values see H&H pg. 72)
MR_design = 1.8# Mixture ratio by weight (ox/fuel)
cstar_eff = 1 #Need to actually implement this


## ADVANCED INPUT PARAMETERS ##
# For fine tuning the nozzle geometry. See Sutton pg. 80 for helpful graphic
adv_data = {"solve_bell": True, # Solve with a bell nozle instead of a conical one
            "div_ang": 15, # Conical divergence half angle [Deg]
            "con_ang": 35, # Convergence half angle [Deg]
            "rad_rat": 0.7, # rad_rat: R/Rmax (BETWEEN 0 & 1)
            "lead_in_factor": 1.5, # Ratio of throat inlet radius of curvature to throat radius
            "lead_out_factor": 0.7, # Ratio of throat outlet radius of curvature to throat radius (has minimal impact)
            "theta_i": 35, # Angle leaving throat [deg] (Should be between 20 and 50)
            "percent_of_conical": 80} # Percent length compared to conical alternative (Should be ~80%)

engine = Engine(thrust, P_c, P_e, con_rat, L_star = L_star, MR = MR_design, adv_data = adv_data, cstar_eff = cstar_eff) # Create engine object
engine.design_engine() # Run engine design process (source code in contourDesigner)
comb_chamber_sizing(engine, axial_vel)

