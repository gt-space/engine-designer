# MAIN.PY – Master file for running engine design scripts
# Calls design.py and gathers results for printing and plotting
import sys
# Uncomment and modify the following sys.path.insert line with the path of your file if "ModuleNotFoundError: No module named '...'" exception is given:
sys.path.insert(0, '/Users/saakethramramoju/Desktop/engineDesigner_v5.0')
import numpy as np
import matplotlib.pyplot as plt
from contourDesigner.design import Engine
import warnings
from heatsinkDesigner.heatsink import Heatsink
from regenDesigner.bartz import bartz
#from regenDesigner.bartz import bartz_x
warnings.filterwarnings("ignore", message="Support for FigureCanvases without a required_interactive_framework attribute was deprecated")
np.set_printoptions(linewidth=1000000)
## GENERAL INPUT PARAMETERS ##
# For basic engine configuration and adjuqstment

thrust = 2500 * 4.4482216153 # Thrust [lbf to N]
P_c = 300 * 0.0689476 # Chamber pressure at injector face [psi to bar]
P_e = 10 / 14.508 # Desired exit presure for expansion ratio [bar] (if in doubt put ambient)
con_rat = 4.3 # Contraction ratio
L_star = 1.02 # Characteristic length [m] (for recommended values see H&H pg. 72)
MR = 2 # Mixture Ratio
cstar_eff = 1.0 # C-star Effeciency
numPTS = 100 # Number of axial stations along contour
fuel = 'JetA' # Fuel type (https://rocketcea.readthedocs.io/en/latest/propellants.html)
ox = 'LOX' # Oxidizer type (https://rocketcea.readthedocs.io/en/latest/propellants.html)
preset_chamber_ID = 5.75 # Chamber inner diameter [in]

'''
# Subscale
thrust = 800 * 4.4482216153 # Thrust [lbf to N]
P_c = 261 * 0.0689476 # Chamber pressure at injector face [psi to bar]
P_e = 18 / 14.508 # Desired exit presure for expansion ratio [bar] (if in doubt put ambient)
con_rat = 4.5 # Contraction ratio
L_star = 1.1 # Characteristic length [m] (for recommended values see H&H pg. 72)
MR = 2 # Mixture Ratio
cstar_eff = 1.0 # C-star Effeciency
numPTS = 200 # Number of axial stations along contour
fuel = 'JetA' # Fuel type (https://rocketcea.readthedocs.io/en/latest/propellants.html)
ox = 'LOX' # Oxidizer type (https://rocketcea.readthedocs.io/en/latest/propellants.html)
preset_chamber_ID = 3.68 # Chamber inner diameter [in]
'''

## For Heatsink Design ##
thickness = 0.85 # Chamber wall thickness [in]
hotfire_time = 10 # Duration of hotfire [s]
graphite_OD = 0 # Outer diameter of graphite insert [in]
graphite_start_index = 126 # The axial station index at which the gprahite insert starts 
graphite_end_index = 151 # The axial station index at which the gprahite insert ends
dt = 1 # Time sub-interval for transient FDM analysis. Make this smaller for finer results [s]
radial_subdivisions = 100 # Number of radial subdivisons along wall thickness at each axial station to perform FDM analysis. Make this larger for finer results 
analysis_index = 1 # Axial station to analyze for 2D plots 

## ADVANCED INPUT PARAMETERS ##
# For fine tuning the nozzle geometry. See Sutton pg. 80 for helpful graphic
adv_data = {"solve_bell": True, # Solve with a bell nozzle instead of a conical one
            "div_ang": 15, # Conical divergence half angle [Deg]
            "con_ang": 35, # Convergence half angle [Deg]
            "rad_rat": 0.7, # rad_rat: R/Rmax (BETWEEN 0 & 1)
            "lead_in_factor": 1.5, # Ratio of throat inlet radius of curvature to throat radius
            "lead_out_factor": 0.7, # Ratio of throat outlet radius of curvature to throat radius (has minimal impact)
            "theta_i": 35, # Angle leaving throat [deg] (Should be between 20 and 50)
            "percent_of_conical": 80} # Percent length compared to conical alternative (Should be ~80%)

## CALL DESIGN SCRIPT ##
preset_chamber_ID = preset_chamber_ID / 39.37 # inches to meters
engine = Engine(thrust, P_c, P_e, con_rat, L_star = L_star, MR = MR, adv_data = adv_data, cstar_eff = cstar_eff, numPTS = numPTS, fuel = fuel, ox = ox, preset_chamber_ID=preset_chamber_ID) # Create engine object given params
engine.design_engine() # Call the design function


# Export contour data to csv
def create_csv(engine):
    matrix = engine.engineContour# Radius and distance from injector matrix
    np.savetxt("contour.csv", matrix, delimiter=",")

create_csv(engine) # Call csv function

print(" +++ ENGINE DESIGN RESULTS +++")
print("Exit Velocity [m/s]: " + str(engine.V_exit))
print("Mass Flow Rate [kg/s]: " + str(engine.mDot_tot))
print("Fuel Mass Flow Rate [kg/s]: " + str(engine.mDot_f))
print("LOX Mass Flow Rate [kg/s]: " + str(engine.mDot_o))
print("Throat Area [m^2]: " + str(engine.A_t))
print("C-star [m/s]: " + str(engine.C_star))
print("Exit angle [deg]: " + str(engine.theta_e))
print("Chamber Radius [in]: " + str(engine.engineProps[0, 0]/.0254))
print("Exit Radius [in]: " + str(engine.engineProps[-1, 0]/.0254))
print("Expansion Ratio: " + str(np.pi * engine.engineProps[-1, 0]**2 / engine.A_t))
print(" ")
plt.figure()
plt.plot(engine.engineProps[:, 1] / 0.0254, engine.engineProps[:, 0] / 0.0254)
plt.xlabel('Distance from Injector [in]', fontsize=16)
plt.ylabel('Radius [in]', fontsize=16)
plt.title("Engine Contour")
plt.xlim(0, engine.engineContour[-1, 1] / .0254)
plt.ylim(0, engine.engineContour[-1, 1] / .0254)
np.savetxt("contour.csv", engine.engineProps[:, 0:2]*39.37, delimiter=',')


heat = Heatsink(engine, thickness = thickness, hotfire_time = hotfire_time, chamber_inner_diameter =  preset_chamber_ID*39.37, graphite_OD = graphite_OD, graphite_start_index = graphite_start_index, graphite_end_index = graphite_end_index, dt = dt, radial_subdivisions = radial_subdivisions)
# Steel and Graphite properties can be modified in heatsink.py

#print(heat.graphite_test()) # determines which axial indices to start and end graphite insert
temps, hg_list = heat.transient_solution() # finds full transient thermal solution at all time and axial indices. All plotting functions automatically do this.
#print(temps) # K, print full wall temperature history (3d array)
#print(hg_list[-1, 1]) # W/m^2-K, print convective coefficient at a given [time index, axial station]
#print(heat.temps[hotfire_time, :, :]) # heat.temps contains full temp history (self.temps[time_index, axial_index, radial_index])
#heat.plot_transient_3d() # full 3d transient plot at each time index
#heat.plot_wall_temp_gradient_at_station(hotfire_time, analysis_index)
#heat.plot_inner_wall_temp_at_time(hotfire_time)
#heat.plot_inner_wall_temp_at_station(analysis_index)
#print(heat.mat) # visual representaion of matrial distribution throughout engine



plt.show()

#np.savetxt("CEA_props_TR7000_PC300_CR3-75_MR2_PE8_LS1-05.csv", engine.engineProps, delimiter=",")
#np.savetxt("Temperature_History.csv", temps, delimiter=",") EDIT THIS