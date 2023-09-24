import numpy as np
import matplotlib.pyplot as plt
from contourDesigner.design import Engine




def mdot_fuel(P_0):
    Cd = 0.76
    A = 2.056 * 10 ** -5
    rho = 800
    P_1 = 25.2 * 100000

    mdot = Cd * A * np.sqrt(2 * rho * np.abs((P_1 - P_0)))
    if P_1 < P_0:
        mdot = mdot * -1
    return mdot

def mdot_ox(P_0):
    Cd = 0.668
    A = 5.255 * 10**-5
    rho = 1141
    P_1 = 21.24 * 100000

    mdot = Cd * A * np.sqrt(2 * rho * np.abs((P_1 - P_0)))
    if P_1 < P_0:
        mdot = mdot * -1
    return mdot

def mdot_nozzle(P_0, engine):
    A_star = engine.A_t

    gam = engine.engineProps[engine.throatInd+100, 15]
    R_sp = 8.315 / engine.engineProps[engine.throatInd + 100, 13] * 1000

    T_0 = engine.engineProps[0, 9]

    return A_star * np.sqrt(gam / R_sp * (2/(gam+1))**((gam+1)/(gam-1))) * P_0 / np.sqrt(T_0)


print(mdot_fuel(1800000) + mdot_ox(1800000))

## GENERAL INPUT PARAMETERS ##
# For basic engine configuration and adjustment
thrust = 786 * 4.44822 # Thrust [lbf to N]
P_c = 261 * 0.0689476 # Chamber pressure at injector face [psi to bar]
P_e = 8 / 14.5038 # Desired exit presure for expansion ratio [bar] (if in doubt put ambient)
con_rat = 5 # Contraction ratio
L_star = 1.05 # Characteristic length [m] (for recommended values see H&H pg. 72)
MR = 1.8# Mixture ratio by weight (ox/fuel)
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

## CALL ENGINE DESIGN SCRIPTS ##
engine = Engine(thrust, P_c, P_e, con_rat, L_star = L_star, MR = MR, adv_data = adv_data, cstar_eff = cstar_eff) # Create engine object
engine.design_engine() # Run engine design process (source code in contourDesigner)

print(mdot_nozzle(1800000, engine))

P_0 = 0
h = 0.000001
time_arr = np.arange(0, 0.1, h)
P_0_arr = np.zeros(len(time_arr) + 1)
R_sp = 8314 / engine.engineProps[0, 13]
T_0 = engine.engineProps[0, 9]
V_c = engine.A_t * L_star
for i in range(len(time_arr)):

    dP0_dt = R_sp * T_0 / V_c * (mdot_fuel(P_0_arr[i]) + mdot_ox(P_0_arr[i]) - mdot_nozzle(P_0_arr[i], engine))
    P_0_arr[i+1] = P_0_arr[i] + dP0_dt * h



plt.figure()
plt.plot(time_arr, P_0_arr[1:]/ 6894.76)
plt.xlabel("Time (sec)")
plt.ylabel("Chamber Pressure (psi)")
plt.title("Transient Startup Pressure")
plt.show()