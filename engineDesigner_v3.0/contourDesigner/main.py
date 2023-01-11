# MAIN.PY – Master file for running engine design scripts
# Calls design.py and gathers results for printing and plotting

from design import Engine
import matplotlib.pyplot as plt
import numpy as np

## GENERAL INPUT PARAMETERS ##
# For basic engine configuration and adjustment
thrust = 4000 * 4.44822 # Thrust [lbf to N]
P_c = 250 * 0.0689476 # Chamber pressure at injector face [psi to bar]
P_e = 1 # Desired exit presure for expansion ratio [bar] (if in doubt put ambient)
con_rat = 5.5 # Contraction ratio
L_star = 1.2 # Characteristic length [m] (for recommended values see H&H pg. 72)
MR = 1.8 # Mixture ratio by weight (ox/fuel)

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
engine.isentropic_comparison() # Compare with isentropic flow relations (optional, prints to console)


## DISPLAY RESULTS ##

print(" +++ RESUTLS +++")
print("Exit Velocity [m/s]: " + str(engine.V_exit))
print("Mass Flow Rate [kg/s]: " + str(engine.mDot_tot))
print("Fuel Mass Flow Rate [kg/s]: " + str(engine.mDot_f))
print("LOX Mass Flow Rate [kg/s]: " + str(engine.mDot_o))
print("Throat Area [m^2]: " + str(engine.A_t))
print("C-star [m/s]: " + str(engine.C_star))
# print("Thrust [lbf]: " + str(engine.thrust * 0.224809)) # Only relevant if using alternative method
print("Exit Pressure [bar]: " + str(engine.engineProps[-1, 8]))
print("Exit angle [deg]: " + str(engine.theta_e))


















# Plot a specific value (change the 0 to the index of the property you want to see):
# Indecies can be found in design.py header
plt.plot(engine.engineProps[:, 1], engine.engineProps[:, 0])
plt.xlabel('Distance from Injector [m]', fontsize=16)
plt.ylabel('Radius [m]', fontsize=16)
# plt.ylim([0,0.1]) # Restrict plot range
plt.show()

# Export contour data to csv
def create_csv(engine):
    matrix = engine.engineProps# Radius and distance from injector matrix
    np.savetxt("CEA_props_6000_400.csv", matrix, delimiter=",")

create_csv(engine) # Call csv function
