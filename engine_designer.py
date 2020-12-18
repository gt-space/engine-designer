from engine_designer.dataCollectionScript import Engine
import matplotlib.pyplot as plt
import numpy as np

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
# vis, cond, pran, condfz, pranfz]

thrust = 17792.89 # thrust (N)
C_p = 48.2633 # chamber pressure (bar)
conrat = 6
LStar = 1.2
# Run it!
engine = Engine(thrust, C_p, conrat, LStar = LStar)
engine.design_engine()
print(engine.engineProps)

print(engine.A_t)
print("C-star:", engine.C_star)

# See a specific value (change the 0 to the index of the property you want to see):
plt.plot(engine.engineProps[:, 1], engine.engineProps[:,0])
plt.xlabel('Distance from Injector (m)', fontsize=16)
plt.ylabel('Radius (m)', fontsize=16)
# plt.ylim([0,0.1])
plt.show()

def create_csv(engine):
    matrix = engine.engineProps[:, 0:2] # Radius and distance from injector matrix
    np.savetxt("contour.csv", matrix, delimiter=",")

# create_csv(engine)
