import numpy as np
from matplotlib import pyplot as plt
from design import Engine
from ablative import Ablative
import warnings
warnings.filterwarnings("ignore", message="Support for FigureCanvases without a required_interactive_framework attribute was deprecated")
plt.close('all')

# ======================= GENERAL INPUT PARAMETERS  ============================
# Subscale
# For basic engine configuration and adjustment
thrust = 800*4.44822 # Thrust [N]
P_c = 18 # Chamber pressure at injector face [psi to bar]
P_e = 12 / 14.508 # Desired exit presure for expansion ratio [bar] (if in doubt put ambient)
#con_rat = 8.0 * D_t^-0.6 + 1.25 where D_t is in cm
con_rat = [] #  Contraction ratio sized using correlation in code
L_star = 1.1 # Characteristic length [m] (for recommended values see H&H pg. 72)
MR = 1.8 # Mixture Ratio

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
engine = Engine(thrust, P_c, P_e, con_rat, L_star = L_star, MR = MR, adv_data = adv_data) # Create engine object given params
engine.design_engine() # Call the design function
engine.isentropic_comparison()
ablative = Ablative(engine, run_time=40, T_wall_i=3220, split=0)

# =======================  Engine Design Results  ============================
print(" +++ ABLATIVE ENGINE DESIGN RESULTS +++")
print("Pre-ablation Exit Velocity [m/s]: " + str(engine.V_exit))
print("Mass Flow Rate [kg/s]: " + str(engine.mDot_tot))
print("Fuel Mass Flow Rate [kg/s]: " + str(engine.mDot_f))
print("LOX Mass Flow Rate [kg/s]: " + str(engine.mDot_o))
print("Pre-ablation Throat Area [in^2]: " + str(engine.A_t/((.0254)**2)))
print("C-star [m/s]: " + str(engine.C_star))
print("Exit angle [deg]: " + str(engine.theta_e))
print("Pre-ablation Chamber Radius [in]: " + str(engine.engineProps[0, 0]/.0254))
print("Pre-ablation Exit Radius [in]: " + str(engine.engineProps[-1, 0]/.0254))
print("Pre-ablation Expansion Ratio: " + str(np.pi * engine.engineProps[-1, 0]**2 / engine.A_t))
print("Pre-ablation Contraction Ratio: " + str(np.pi * engine.engineProps[0, 0]**2 / engine.A_t))
print("Engine Length [in]: " + str(engine.engineProps[-1, 1]/.0254))

#print("Chamber Length [in]:"+ str(engine.engineProps[engine.throat_end_ind, 1]/.0254)) This isn't working rn, using plot

print(" ")

# ======================= Engine Contour Plot  ============================

plt.figure()
plt.plot(engine.engineProps[:, 1] / 0.0254, engine.engineProps[:, 0] / 0.0254)
plt.xlabel('Distance from Injector [in]', fontsize=16)
plt.ylabel('Radius [in]', fontsize=16)
plt.title("Engine Contour")
plt.xlim(0, engine.engineContour[-1, 1] / .0254)
plt.ylim(0, engine.engineContour[-1, 1] / .0254)
plt.show()
# Export Contour to CAD
np.savetxt('CONTOUR.txt', engine.engineContour, delimiter=' ', fmt='%s')
# ======================= Engine Mass Analysis  ============================
# Initial Inputs
t_w = 1*.0254
t_I = 0*.0254
t_OW =0.125*.0254
mat_w= Ablative.SilicaEpoxy()
mat_I=ablative.SilicaEpoxy()
mat_OW=ablative.CarbonEpoxy()
mass = ablative.estim_liner_mass(engine,t_w,t_I,t_OW,mat_w,mat_I,mat_OW)
print("Engine mass (no injector) [lbm]: " + str(mass))

# ======================= Throat Insert  ============================
t_th = ablative.throat_insert_size(ablative.Graphite(),3)
print("Throat Insert Thickness [in]: " + str(t_th/.0254))

# Tolerance
'''# This is used to rationalize ten-thou tolerance on graphite throat insert
tol=.01*.0254 #m
thrust_new= ablative.thrust_change_throat(tol,0, 1)[0]
print("===Tolerance Analysis===")
print("Design thrust (N): " + str(thrust))
print("New thrust (N): " + str(thrust_new) )
print("% Change thrust: " + str((thrust_new-thrust)/thrust*100))'''

# ======================= Chamber Pressure Drop Modeling  ============================
# This is a first-order model of how chamber pressure dropping could affect thrust
# This is incorrect; I don't expect the chamber pressure drop to be this severe, but it is useful to see
# how much thrust will change as a result of chamber pressure changing
'''
chamber_ablation_rate = .0007*.0254  # guess from RDL numbers (m/s)
throat_insert_ablation_rate = .0006 *.0254 # Sutton

p_c_f = ablative.p_c_drop(chamber_ablation_rate,40,showplot=True)

burnTime = 40 #seconds


print("")
print("===Ablative Chamber Pressure Drop Analysis===")
print("Final Chamber Pressure (psi) " + str(p_c_f/6894.76))

ablative.thrust_curve(throat_insert_ablation_rate,chamber_ablation_rate, burnTime,debugMach=False,showplot=True)
'''

'''# ======================= Overwrap Structural Analysis  ============================
#Material Selection
CF = ablative.CarbonEpoxy()
liner_mat = ablative.CarbonEpoxy()

# Geometry Choices
t_OW = .125*.0254
t_liner = .5*.0254
T_amb = 300 #K

(t_A, t_OW, q_conv, h_g, T_aw, T_I, T_OW)=ablative.size_thickness(t_OW)
DT_liner = T_aw-T_I
DT_OW = T_I-T_amb
(stress_OW,FOS_OW)= ablative.chamber_stress_analysis(CF,t_OW,DT_OW)
(stress_liner,FOS_liner)= ablative.chamber_stress_analysis(liner_mat,t_liner,DT_liner)

print("")
print("===Chamber Structural Analysis===")
#print("Liner Stress: " + str(stress_liner))
#print("Liner FOS: " + str(FOS_liner))
# Liner stresses are not relevant because the DT should never get to steady state; ablation should start before then
print("OW Stress: " + str(stress_OW))
print("Temp Gradient across Outer Wall (K): " + str(DT_OW))
print("OW FOS: " + str(FOS_OW))
plt.show()
'''