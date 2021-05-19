# MAIN.PY – Master file for running full engine analysis
# Effectively the same as main.py in contourDesigner, but with regen functionality

from contourDesigner.design import Engine
from regenDesigner.regen import RegenJacket
import matplotlib.pyplot as plt
import numpy as np

## GENERAL INPUT PARAMETERS ##
# For basic engine configuration and adjustment
thrust = 4000 * 4.44822 # Thrust [lbf to N]
P_c = 250 * 0.0689476 # Chamber pressure at injector face [psi to bar]
P_e = 0.5 # Desired exit presure for expansion ratio [bar] (if in doubt put ambient)
con_rat = 5 # Contraction ratio
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

## REGEN PARAMETERS ##
channel_height = 0.002 # Height of channel (depth of endmill cut) [m]
wall_thickness = 0.002 # Thicness of liner (distance from inner wall to channel bottom) [m]
min_fin_width = 0.001 # Fin width at the throat [m]
channel_width = 0.002 # Width of channel [m]
starting_index = 45 # Index to begin analysis for non-full regen designs. Where the fuel enters


## CALL ENGINE DESIGN SCRIPTS ##
engine = Engine(thrust, P_c, P_e, con_rat, L_star = L_star, MR = MR, adv_data = adv_data) # Create engine object
engine.design_engine() # Run engine design process (source code in contourDesigner)

## DISPLAY RESULTS ##
print(" +++ ENGINE DESIGN RESUTLS +++")
print("Exit Velocity [m/s]: " + str(engine.V_exit))
print("Mass Flow Rate [kg/s]: " + str(engine.mDot_tot))
print("Fuel Mass Flow Rate [kg/s]: " + str(engine.mDot_f))
print("LOX Mass Flow Rate [kg/s]: " + str(engine.mDot_o))
print("Throat Area [m^2]: " + str(engine.A_t))
print("C-star [m/s]: " + str(engine.C_star))
# print("Thrust [lbf]: " + str(engine.thrust * 0.224809)) # Only relevant if using alternative method
print("exit pressure [bar]: " + str(engine.engineProps[-1, 8]))
print("Exit angle [deg]: " + str(engine.theta_e))
print(" ")

# Plot a specific value (change the 0 to the index of the property you want to see):
# Indecies can be found in design.py header
plt.plot(engine.engineProps[:, 1], engine.engineProps[:, 0])
plt.xlabel('Distance from Injector [m]', fontsize=16)
plt.ylabel('Radius [m]', fontsize=16)
plt.show()

print("Running Regen Analysis...")

## CALL REGEN SCRIPTS ##
jacket = RegenJacket(engine, channel_height, wall_thickness, min_fin_width, channel_width, start_ind = starting_index) # Create jacket object
(wall_temps, coolant_temps, pressures, num_channels) = jacket.simulate_regen() # Run regen simulation
# Gets wall temperatures, coolant temperatures, pressures, and channel number

## DISPLAY RESULTS ##
print(" +++ REGEN SIM RESULTS +++")
print("Number of Channels: " + str(num_channels))
print("Max wall temp [K]: " + str(max(wall_temps)[0]))
print("Pressure Loss [psi]: " + str((pressures[-1]-pressures[0])[0]*14.5038))
print("Coolant Outlet Temperature [K]: " + str(coolant_temps[0][0]))

## PLOT REGEN RESULTS ##
# Change wall_temps to the vector of your chosing
plt.plot(engine.engineProps[:-(starting_index + 1):,1], wall_temps)
plt.xlabel('Distance from Injector (m)', fontsize=16)
plt.ylabel('Wall Temperature (K)', fontsize=16)
plt.show()
