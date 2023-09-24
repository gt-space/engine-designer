# ENGINE_TRADES.PY – spinoff from main_KeraLOX.py for analyzing combinations of engine configurations
# Effectively the same as main_KeraLOX.py in contourDesigner, but with regen functionality

from contourDesigner.design import Engine
from regenDesigner.regen import RegenJacket
import matplotlib.pyplot as plt

import numpy as np
from bartz import bartz
from chan_contour import get_chan_contour


## GENERAL INPUT PARAMETERS ##
# For basic engine configuration and adjustment
thrust = 1000 * np.array([22, 24, 26, 28, 30, 32])  #* 4.44822 # Thrust [lbf to N]
P_c = [18, 21, 24, 27, 30, 33, 36] # 0.0689476 # Chamber pressure at injector face [psi to bar]
P_e = 8 / 14.5038 # Desired exit presure for expansion ratio [bar] (if in doubt put ambient)
con_rat = 4 # Contraction ratio
L_star = 1.05 # Characteristic length [m] (for recommended values see H&H pg. 72)
MR_design = 2.3 #Mixture ratio by weight (ox/fuel)
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

## REGEN PARAMETERS ##
wall_thickness = 0.001 # Thickness of liner (distance from inner wall to channel bottom) [m]
min_fin_width = 0.00075 # Fin width at the throat [m]
min_channel_width = 0.0005 # Minimum allowable channel width [m] (determined by manufacturing capability)
var_channel_width = False # Turn variable channel width on or off (e.g., using an end mill vs a slitting saw)
starting_index = 0 # Index to begin analysis for non-full regen designs. Where the fuel enters

min_channel_Re = 4000 # Minimum Reynolds number in channels
max_fin_AR = 3 # Maximum fin aspect ratio (height/width)
max_channel_h = 0.01905 # (m) Max channel height
#max_channel_h = 0.004

outlet_temp = 550 # (K) Max coolant outlet temp
uts_FOS_min = 2 # Minimum required ultimate factor of safety
max_wall_temp = 800 # (K)

#LOX_Mdot = np.zeros([len(P_c), len(thrust)])
#Fuel_Mdot = np.zeros([len(P_c), len(thrust)])

exit_diam = np.zeros([len(P_c), len(thrust)])

for i in np.arange(len(P_c)):
    for j in np.arange(len(thrust)):
        ## CALL ENGINE DESIGN SCRIPTS ##
        engine = Engine(thrust[j], P_c[i], P_e, con_rat, L_star = L_star, MR = MR_design, adv_data = adv_data, cstar_eff = cstar_eff) # Create engine object
        engine.design_engine() # Run engine design process (source code in contourDesigner)

        #LOX_Mdot[i, j] = engine.mDot_o
        #Fuel_Mdot[i, j] = engine.mDot_f

        exit_diam[i, j] = engine.engineProps[-1, 0]


#np.savetxt('CEA_781_260_8.csv', engine.engineProps, delimiter=",")

## DISPLAY RESULTS ##
print(" +++ ENGINE DESIGN RESULTS +++")
print("Exit Velocity [m/s]: " + str(engine.V_exit))
print("Mass Flow Rate [kg/s]: " + str(engine.mDot_tot))
print("Fuel Mass Flow Rate [kg/s]: " + str(engine.mDot_f))
print("LOX Mass Flow Rate [kg/s]: " + str(engine.mDot_o))
print("Throat Area [m^2]: " + str(engine.A_t))
print("C-star [m/s]: " + str(engine.C_star))
print("Engine ISP [s]: " + str(engine.V_exit / 9.8104))
print("Exit angle [deg]: " + str(engine.theta_e))
print("Chamber Radius [m]: " + str(engine.engineProps[0, 0]))
print("Exit Radius [m]: " + str(engine.engineProps[-1, 0]))
print("Expansion Ratio: " + str(np.pi * engine.engineProps[-1, 0]**2 / engine.A_t))
print(" ")