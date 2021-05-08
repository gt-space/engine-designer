# MAIN.PY – Master file for running engine design scripts
# Calls design.py and gathers results for printing and plotting

from design import Engine
from alt_design import AltEngine
import matplotlib.pyplot as plt
import numpy as np

## DEFINE INPUT PARAMETERS ##
thrust = 4500 * 4.44822 # Thrust [lbf to N]
P_c = 500 * 0.0689476 # Chamber pressure at injector face [psi to bar]
con_rat = 6 # Contraction ratio
L_star = 1.2 # Characteristic length [m]
MR = 1.8 # Mixture ratio by weight (ox/fuel)

# Params for alternative solution method (uses throat area instead of thrust)
D_exit = 10 * 0.0254 # Based on guess on vehicle diameter [in to m]
alt = False # Use alternative solution method where thrust is solved for

## CALL DESIGN SCRIPT ##
if alt == False:
    engine = Engine(thrust, P_c, con_rat, L_star = L_star, MR = MR) # Create engine object given params
    engine.design_engine() # Call the design function
else:
    engine = AltEngine(D_exit, P_c, con_rat, L_star = L_star, MR = MR) # Create engine object given params
    engine.design_engine() # Call the design function

## COMPARE WITH ISENTROPIC RELATIONS ##

engine.isentropic_comparison()

## DISPLAY RESULTS ##
print("Exit Velocity [m/s]: " + str(engine.V_exit))
print("Mass Flow Rate [kg/s]: ", str(engine.mDot_tot))
print("Fuel Mass Flow Rate [kg/s]: " + str(engine.mDot_f))
print("LOX Mass Flow Rate [kg/s]: " + str(engine.mDot_o))
print("Throat Area [m^2]: " + str(engine.A_t))
print("C-star [m/s]: ", str(engine.C_star))
print("Thrust [lbf]: ", str(engine.thrust * 0.224809))
print("exit pressure [bar]: ", str(engine.engineProps[-1, 8]))

# Plot a specific value (change the 0 to the index of the property you want to see):
# Indecies can be found in design.py header
plt.plot(engine.engineProps[:, 1], engine.engineProps[:, 2])
plt.xlabel('Distance from Injector [m]', fontsize=16)
plt.ylabel('Radius [m]', fontsize=16)
# plt.ylim([0,0.1]) # Restrict plot range
plt.show()

# Export contour data to csv
def create_csv(engine):
    matrix = engine.engineProps[:, 0:2] # Radius and distance from injector matrix
    np.savetxt("contour.csv", matrix, delimiter=",")

# create_csv(engine) # Call csv function
