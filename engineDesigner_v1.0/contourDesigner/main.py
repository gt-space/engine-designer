# USAGE:
# design_engine(thrust_nom, P_inj, conRat, MR, divAng)
# thrust_nom: Design thrust (N)
# P_inj: Injector face pressure (Bar)
# conRat: Enigne contraciton ratio
# MR: ox/fuel (default = 1.8)
# divAng: Divergence half angle (Deg) (defualt = 15)

# OUTPUTS:
# engineContour is a numpy array with the following column format:
# [R, Z]
# engineProps is a numpy array with the following column format:
# [R, Z, pip, aeat, mach, cf, ivac, isp, p, t, rho, h, u, mw, cp, gam, son,
# vis, cond, pran, cpfz, condfz, pranfz]

from engine_designer.dataCollectionScript import Engine
import matplotlib.pyplot as plt
import numpy as np

# Define input parameters
thrust = 3000 * 4.44822 # Thrust [lbf to N]
P_c = 500 * 0.0689476 # Chamber Pressure [psi to bar]
conrat = 6 # Contraction Ratio
LStar = 1.2 # Characteristic Length [m]
MR = 1.8 # Mixture ratio by weight

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
