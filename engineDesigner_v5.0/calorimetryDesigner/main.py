# MAIN.PY – Cooling channel optimization with real engine geometry

import sys
sys.path.insert(0, '/Users/kieranyarberry/Desktop/engineDesigner_v5.0')
import numpy as np
from contourDesigner.design import Engine
from calChamber.calChamber2 import CoolingChannelDesigner, MaterialProperties
import matplotlib.pyplot as plt
from calChamber.calChamber2 import CoolingChannelDesigner
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
from plotting import plot_engine_3d, plot_segment_profiles

# ENGINE DESIGN PARAMETERS
thrust = 3200 * 4.4482216153  # Thrust [lbf to N]
P_c = 400 * 0.0689476         # Chamber pressure at injector face [psi to bar]
P_e = 10 / 14.508             # Desired exit pressure [bar]
con_rat = 4.3                 # Contraction ratio
L_star = 1.02                 # Characteristic length [m]
MR = 2.23                        # Mixture Ratio
cstar_eff = 0.9               # C-star Efficiency
numPTS = 200                  # Number of axial stations
fuel = 'JetA'
ox = 'LOX'
preset_chamber_ID = 5.75 / 39.37  # in → m

# ADVANCED GEOMETRY OPTIONS 
adv_data = {
    "solve_bell": True,
    "div_ang": 15,
    "con_ang": 35,
    "rad_rat": 0.7,
    "lead_in_factor": 1.5,
    "lead_out_factor": 0.7,
    "theta_i": 35,
    "percent_of_conical": 80
}

# ENGINE INIT AND DESIGN 
engine = Engine(
    thrust, P_c, P_e, con_rat,
    L_star=L_star, MR=MR,
    adv_data=adv_data, cstar_eff=cstar_eff,
    numPTS=numPTS, fuel=fuel, ox=ox,
    preset_chamber_ID=preset_chamber_ID
)
engine.design_engine()

# MATERIAL PROPERTIES 
material = MaterialProperties(
    conductivity=40,
    yield_strength=2.48e8,
    ultimate_strength=3.00e8
)

# USER SETTINGS
n_segments = 16
settings = {
    "n_segments": n_segments,
    "n_azimuthal_divisions": 20,
    "P_inlet": 2.758e6,          # Pa
    "T_in": 300,             # K
    "Target_Wall_Temp": 1500,  # K
    "inner_wall_thickness": 0.0015,  # m
    'min_FOS': 1.5,  
    "max_dT": None,  # K, Max allowable (optional)
    "max_dP": None,  # Pa, Max allowable (optional)
    "Superheat": 5, # K, Max allowable wall superheat, 5K to prevent onset of nucleate boiling 

    "geometry_controls": [
        {
            "fix": [],
            "constraints": {
                "min_h": 0.001,
                "max_h": 0.1,
                "min_w": 0.0005,
                "max_w": 0.1,
                "min_fin_thickness": 0.00254,
                "max_aspect_ratio": 7.0
            }
        }
    ] * n_segments,

    #Recommended to sweep through a wide range of parameters, and then narrow search down based on constraints and segment results.

    "sweep_ranges": {
        "h": [0.0063], #np.linspace(0.001, 0.00635, 10),
        "w": [0.0011], #np.linspace(0.001, 0.0015875, 10),
        "N": [4],
        "mdot": np.linspace(0.5, 0.1, 20)
    },

    "weights": {
        "dp": 0.5,
        "dT": 1,
        "mdot": 120,
        "Tw": 100
    }
}

# RUN OPTIMIZATION 
designer = CoolingChannelDesigner(engine, material, settings)
result = designer.optimize()

# OUTPUT 
designer.summarize()
designer.export_results("cooling_design_output.csv")

np.savetxt("contour.csv", engine.engineProps[:, 0:2]*39.37, delimiter=',')
#plot_engine_3d(designer, nx=500, ny=360, az_interp='cubic', axial_temp_interp='linear', log_temp=False, elev=100, azim=0)
plot_segment_profiles(designer, seg_index=12, az_oversample=5, az_interp='cubic', pressure_units='Psi', log_temp=False)