from engine_designer.dataCollectionScript import Engine
import matplotlib.pyplot as plt
import numpy as np

# Define input parameters
thrust = 3000 * 4.44822 # Thrust [lbf to N]
P_c = 500 * 0.0689476 # Chamber pressure at injector face [psi to bar]
conrat = 6 # Contraction ratio
LStar = 1.2 # Characteristic length [m]
MR = 1.8 # Mixture ratio by weight (ox/fuel)

# Run it!
engine = Engine(thrust, P_c, conrat, LStar = LStar) # Create engine object given params
engine.design_engine() # Call the design function

# Print some of the results
print("Exit Velocity [m/s]: " + str(engine.V_exit))
print("Total Mass Flow Rate [kg/s]: " + str(engine.mDot_tot))
print("Throat Area [m^2]: " + str(engine.A_t))
print("C-star [m/s]: ", str(engine.C_star))

# See a specific value (change the 0 to the index of the property you want to see):
plt.plot(engine.engineProps[:, 1], engine.engineProps[:,0])
plt.xlabel('Distance from Injector (m)', fontsize=16)
plt.ylabel('Radius (m)', fontsize=16)
# plt.ylim([0,0.1])
plt.show()

# Export contour data to csv
def create_csv(engine):
    matrix = engine.engineProps[:, 0:2] # Radius and distance from injector matrix
    np.savetxt("contour.csv", matrix, delimiter=",")

# create_csv(engine)
