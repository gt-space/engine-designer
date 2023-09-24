import numpy as np
from matplotlib import pyplot as plt
from design import Engine
from ablative import Ablative
from bartz import bartz
import warnings
warnings.filterwarnings("ignore", message="Support for FigureCanvases without a required_interactive_framework attribute was deprecated")
#plt.close('all')

# ======================= GENERAL INPUT PARAMETERS  ============================
# For basic engine configuration and adjustment
vehicle = "Darcy Space"
#vehicle = "KeraLOX"

if vehicle=="Darcy Space":
    thrust = 1400*4.44822 # Thrust [N]
    P_c = (400)*0.0689476 # Chamber pressure at injector face [bar]
    P_e = (8)*0.0689476 # Desired exit presure for expansion ratio [bar] (if in doubt put ambient)
    #con_rat = 8.0 * D_t^-0.6 + 1.25 where D_t is in cm
    con_rat = 5.5
    L_star = 1.5 # Characteristic length [m] (for recommended values see H&H pg. 72)
    # This came from https://www.halfcatrocketry.com/halfcat
    # as well as Darcy II
    # I cut down L* to fit within the envelope
    MR = 4 # Mixture Ratio
    ox = 'N2O'
    fuel = 'Isopropanol'
    t_A = 0.75 * .0254
    t_I = 0.1 * .0254
    t_OW = .1 * .0254

    FOS_th = 2 # this is for the throat insert
    FOS_OW = 3
    insert_contour = False

    #Envelope is 17 in length and 5.53in diameter OD max right now
if vehicle =="KeraLOX":
    thrust = 800 * 4.44822  # Thrust [N]
    P_c = 18  # Chamber pressure at injector face [psi to bar]
    P_e = 12 / 14.508  # Desired exit presure for expansion ratio [bar] (if in doubt put ambient)
    # con_rat = 8.0 * D_t^-0.6 + 1.25 where D_t is in cm
    con_rat = []  # Contraction ratio sized using correlation in code
    L_star = 1.1  # Characteristic length [m] (for recommended values see H&H pg. 72)
    MR = 1.8  # Mixture Ratio
    ox="LOX"
    fuel="JetA"
    t_A = 0.75 * .0254
    t_I = 0
    t_OW = 0.125 * .0254

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
engine = Engine(thrust, P_c, P_e, con_rat, L_star = L_star, MR = MR, adv_data = adv_data,ox = ox,fuel = fuel, cstar_eff=.9) # Create engine object given params
engine.design_engine() # Call the design function

ablative = Ablative(engine, run_time=40, T_wall_i=3220, split=0)

# =======================  Engine Design Results  ============================
print(" +++ ABLATIVE ENGINE DESIGN RESULTS +++")
print("Pre-ablation Exit Velocity [m/s]: " + str(engine.V_exit))
print("Mass Flow Rate [kg/s]: " + str(engine.mDot_tot))
print("Fuel Mass Flow Rate [kg/s]: " + str(engine.mDot_f))
print("Oxidizer Mass Flow Rate [kg/s]: " + str(engine.mDot_o))
print("Pre-ablation Throat Area [in^2]: " + str(engine.A_t/((.0254)**2)))
print("Pre-ablation Throat Radius [in] " +str(engine.R_t/.0254))
print("C-star [m/s]: " + str(engine.C_star))
print("Exit angle [deg]: " + str(engine.theta_e))
print("Pre-ablation Chamber Radius [in]: " + str(engine.engineProps[0, 0]/.0254))
print("Pre-ablation Exit Radius [in]: " + str(engine.engineProps[-1, 0]/.0254))
print("Pre-ablation Expansion Ratio: " + str(np.pi * engine.engineProps[-1, 0]**2 / engine.A_t))
print("Pre-ablation Contraction Ratio: " + str(np.pi * engine.engineProps[0, 0]**2 / engine.A_t))
print("Engine Length [in]: " + str(engine.engineProps[-1, 1]/.0254))

print("Real Throat Radius [in] " +str(engine.R_t/.0254))
#print("Chamber Length [in]:"+ str(engine.engineProps[engine.throat_end_ind, 1]/.0254)) This isn't working rn, using plot
print("ThroatRad 2 (in): " + str(engine.nozzleContour[engine.throatInd,0]/.0254))
print("ThroatRad 3 (in): " + str(engine.engineContour[engine.throatInd_engprops,0]/.0254))
print("CH barrel end: " + str(engine.engineContour[engine.chBarrel_endInd,1]/.0254))
print("")

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
np.savetxt('DARCYCONTOUR_in.txt', engine.engineContour/.0254, delimiter=' ', fmt='%s')



# ======================= Engine Liner Thicknesses  ============================
overwrap_mat = ablative.CarbonEpoxy()
overwrap = ablative.Overwrap(overwrap_mat,FOS=FOS_OW,DT=45,p=400*6894.76,r=engine.engineProps[0, 0])

if vehicle=="Darcy Space":
    ablator_mat = ablative.CarbonEpoxy()
    #ablator_mat = ablative.SilicaEpoxy()
    insulator_mat = ablative.SilicaEpoxy()
    insulator = ablative.Insulator(insulator_mat)
    #insulator = []

if vehicle == "KeraLOX":
    ablator_mat = ablative.SilicaEpoxy()
    insulator =[]

ablator = ablative.Ablator(ablator_mat,burn_time=40,t_A_design=t_A)
burn_time = 40 #s
ablative.run_liner_sizer(overwrap,insulator,ablator,burn_time,t_A_design=t_A*.0254)

# ======================= Throat Insert  ============================
#t_th = ablative.throat_insert_size(ablative.Graphite(),FOS=FOS_th)

#Tolerance
'''# This is used to rationalize ten-thou tolerance on graphite throat insert
tol=.01*.0254 #m
thrust_new= ablative.thrust_change_throat(tol,0, 1)[0]
print("===Tolerance Analysis===")
print("Design thrust (N): " + str(thrust))
print("New thrust (N): " + str(thrust_new) )
print("% Change thrust: " + str((thrust_new-thrust)/thrust*100))'''

# ======================= Engine Mass Analysis  ============================
# Initial Inputs

CC_rad= engine.engineProps[0, 0]
R_CC = t_A+t_I+t_OW+CC_rad
OD_CC = 2*R_CC #m
#t_th = .5*.0254
#len_th = 3.5*.0254
mat_A= ablative.CarbonEpoxy()
mat_I=ablative.SilicaEpoxy()
mat_OW=ablative.CarbonEpoxy()
insert_mat = ablative.Graphite()

insert= ablative.Insert(ablative, insert_mat,FOS=FOS_th,start_ind=[],end_ind=[],contoured=insert_contour)
t_th_mass=insert.mass*2.20462
t_th=insert.t_w


wall_idx=engine.throatInd_engprops
T_w=  engine.engineProps[wall_idx, 9]
(h_g, q_conv, T_aw) = bartz(engine, T_w, wall_idx)

print("T_w (K): " +str(T_w))
print("h_g_throat (W/m^2K) : " +str(h_g))
print("q_conv_throat (W/m^2) : " +str(q_conv))
print("T_th (K): " +str(T_aw))

stress_perc_th = overwrap.stress_dist()

print("Perc of overwrap stresses are thermal: " +str(stress_perc_th*100)+" %")
#print("Throat Insert Thickness [in]: " + str(t_th/.0254))
#print("Throat Insert Mass [lbm]: " + str(t_th_mass))
#print(" ")
#mass = ablative.estim_engine_mass(t_A,t_I,t_OW,mat_A,mat_I,mat_OW,insert)
#print("Engine mass (no injector and no throat insert) [lbm]: " + str(mass))
#ablative.contour(engine,ablator,t_I,overwrap,insert,contoured=insert_contour).print_insert()


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

(t, t_OW, q_conv, h_g, T_aw, T_I, T_OW)=ablative.size_thickness(t_OW)
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