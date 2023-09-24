# MAIN.PY – Master file for running full engine analysis
# Effectively the same as main_KeraLOX.py in contourDesigner, but with regen functionality

from contourDesigner.design import Engine
from regenDesigner.regen import RegenJacket
import matplotlib.pyplot as plt

import numpy as np
from bartz import bartz
from chan_contour import get_chan_contour


## GENERAL INPUT PARAMETERS ##
# For basic engine configuration and adjustment
thrust = 22859 #* 4.44822 # Thrust [lbf to N]
P_c = 250/14.5038 # 0.0689476 # Chamber pressure at injector face [psi to bar]
P_e = 8 / 14.5038 # Desired exit presure for expansion ratio [bar] (if in doubt put ambient)
con_rat = 4 # Contraction ratio
#con_rat = 8.0 * D_t^-0.6 + 1.25 where D_t is in cm
L_star = 1.05 # Characteristic length [m] (for recommended values see H&H pg. 72)
MR_design = 2# Mixture ratio by weight (ox/fuel)
cstar_eff = 1 #Need to actually implement this
numPTS = 250



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
max_channel_h = 0.006 # (m) Max channel height
#max_channel_h = 0.004

outlet_temp = 550 # (K) Max coolant outlet temp
uts_FOS_min = 1.5 # Minimum required ultimate factor of safety
max_wall_temp = 1000 # (K)

## CALL ENGINE DESIGN SCRIPTS ##
engine = Engine(thrust, P_c, P_e, con_rat, L_star = L_star, MR = MR_design, adv_data = adv_data, cstar_eff = cstar_eff, numPTS = numPTS) # Create engine object
engine.design_engine() # Run engine design process (source code in contourDesigner)

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






# Plot a specific value (change the 0 to the index of the property you want to see):
# Indecies can be found in design.py header
plt.figure()
plt.plot(engine.engineProps[:, 1] / .0254, engine.engineProps[:, 0] / .0254)
plt.xlabel('Distance from Injector [in]')
plt.ylabel('Radius [in]')
plt.title('Engine Contour')
plt.xlim(0, engine.engineContour[-1, 1] / .0254)
plt.ylim(0, engine.engineContour[-1, 1] / .0254)
#plt.show()

plt.figure()
plt.plot(engine.engineProps[:, 1] * 1000, engine.engineProps[:, 0] * 1000)
plt.xlabel('Distance from Injector [mm]')
plt.ylabel('Radius [mm]')
plt.title('Engine Contour')
plt.xlim(0, engine.engineContour[-1, 1] * 1000)
plt.ylim(0, engine.engineContour[-1, 1] * 1000)
#plt.show()

num_channels_arr =  np.arange(60,180,10)
channel_w_res = 0.0001

#Optimize regen design
jacket = RegenJacket(engine, wall_thickness, min_fin_width, min_channel_width = min_channel_width,
                     var_channel_width = var_channel_width,T_ci=300, start_ind = starting_index, min_channel_Re=min_channel_Re,
                     max_fin_AR = max_fin_AR, max_channel_h = max_channel_h, max_wall_temp = max_wall_temp)# Create jacket object
MR_run = jacket.MR_optimization(outlet_temp, uts_FOS_min)
jacket.simulate_regen(MR_run, num_channels_arr, channel_w_res)
num_channels = 120
channel_w = 0.0015
jacket.single_design(MR_run, num_channels, channel_w, plot=True)


## DISPLAY RESULTS ##
print(" +++ REGEN SIM RESULTS +++")
print("Number of Channels: " + str(num_channels))
print("Max wall temp [K]: " + str(max(jacket.wall_temps)[0]))
print("Pressure Loss [psi]: " + str((jacket.pressures[-1]-jacket.pressures[0])[0]*14.5038))
print("Coolant Outlet Temperature [K]: " + str(jacket.coolant_temps[0][0]))


get_chan_contour(jacket)

plt.show()





