import numpy as np
from matplotlib import pyplot as plt
from design import Engine
from ablative import Ablative
import warnings
warnings.filterwarnings("ignore", message="Support for FigureCanvases without a required_interactive_framework attribute was deprecated")



# ======================= GENERAL INPUT PARAMETERS  ============================
# Subscale
# For basic engine configuration and adjustment
plt.close('all')

thrust = 800*4.44822 # Thrust [N]
P_c = 18*14.5038 * 0.0689476 # Chamber pressure at injector face [psi to bar]
P_e = 12 / 14.508 # Desired exit presure for expansion ratio [bar] (if in doubt put ambient)
#con_rat = 8.0 * D_t^-0.6 + 1.25 where D_t is in cm
con_rat = 4.5 #  Contraction ratio sized using correlation in code
L_star = 1.1 # Characteristic length [m] (for recommended values see H&H pg. 72)
MR = 2 # Mixture Ratio

'''# Daedelus

thrust = 19500 # Thrust [N]
P_c = 250 * 0.0689476 # Chamber pressure at injector face [psi to bar]
P_e = 12 / 14.508 # Desired exit presure for expansion ratio [bar] (if in doubt put ambient)
#con_rat = 8.0 * D_t^-0.6 + 1.25 where D_t is in cm
con_rat = 4.5 #  Contraction ratio sized using correlation in code
L_star = 1.1 # Characteristic length [m] (for recommended values see H&H pg. 72)
MR = 2 # Mixture Ratio'''


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

# ======================= Transient Heat Transfer Analysis  ============================
run_time = 10 #secs
T_wall_i = 294

split = 0
ablative = Ablative(engine, run_time, T_wall_i, split)

wall_temp_arr = ablative.chamber_analysis(run_time,liner_thickness=.25*.0254, iter = 3, plot2d=True, plot3d = True)
plt.show()

T_arr_2d =  wall_temp_arr[-1, :]

plt.show()

np.savetxt("CEA_props_TR7000_PC300_CR3-75_MR2_PE8_LS1-05.csv", engine.engineProps, delimiter=",")
