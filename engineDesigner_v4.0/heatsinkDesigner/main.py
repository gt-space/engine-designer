import numpy as np
from matplotlib import pyplot as plt
from design import Engine
from heatsink import HeatSink

## GENERAL INPUT PARAMETERS ##
# For basic engine configuration and adjustment
thrust = 22859 # Thrust [lbf to N]
P_c = 250 * 0.0689476 # Chamber pressure at injector face [psi to bar]
P_e = 8 / 14.508 # Desired exit presure for expansion ratio [bar] (if in doubt put ambient)
con_rat = 3.75 # Contraction ratio
L_star = 1.05 # Characteristic length [m] (for recommended values see H&H pg. 72)
MR = 2.0 # Mixture Ratio

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

## CALL DESIGN SCRIPT ##
engine = Engine(thrust, P_c, P_e, con_rat, L_star = L_star, MR = MR, adv_data = adv_data) # Create engine object given params
engine.design_engine() # Call the design function

print(" +++ ENGINE DESIGN RESULTS +++")
print("Exit Velocity [m/s]: " + str(engine.V_exit))
print("Mass Flow Rate [kg/s]: " + str(engine.mDot_tot))
print("Fuel Mass Flow Rate [kg/s]: " + str(engine.mDot_f))
print("LOX Mass Flow Rate [kg/s]: " + str(engine.mDot_o))
print("Throat Area [m^2]: " + str(engine.A_t))
print("C-star [m/s]: " + str(engine.C_star))
print("Exit angle [deg]: " + str(engine.theta_e))
print("Chamber Radius [m]: " + str(engine.engineProps[0, 0]))
print("Exit Radius [m]: " + str(engine.engineProps[-1, 0]))
print("Expansion Ratio: " + str(np.pi * engine.engineProps[-1, 0]**2 / engine.A_t))
print(" ")

plt.plot(engine.engineProps[:, 1] / 0.0254, engine.engineProps[:, 0] / 0.0254)
plt.xlabel('Distance from Injector [in]', fontsize=16)
plt.ylabel('Radius [in]', fontsize=16)
plt.title("Engine Contour")
plt.xlim(0, engine.engineContour[-1, 1] / .0254)
plt.ylim(0, engine.engineContour[-1, 1] / .0254)

#Hot Fire parameters
run_time = 40 #secs
T_wall_i = 294

split = 0
heatsink = HeatSink(engine, run_time, T_wall_i, split)

plt.figure()
wall_temp = heatsink.chamber_analysis(run_time, iter = 5, plot2d=True, plot3d = True)

T_arr_2d =  wall_temp[-1, :]

plt.show()

np.savetxt("CEA_props_TR7000_PC300_CR3-75_MR2_PE8_LS1-05.csv", engine.engineProps, delimiter=",")