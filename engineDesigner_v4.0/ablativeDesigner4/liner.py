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
P_c = 18*14.5038 * 0.0689476 # Chamber pressure at injector face [psi to bar]
P_e = 12 / 14.508 # Desired exit presure for expansion ratio [bar] (if in doubt put ambient)
#con_rat = 8.0 * D_t^-0.6 + 1.25 where D_t is in cm
con_rat = 4.5 #  Contraction ratio sized using correlation in code
L_star = 1.1 # Characteristic length [m] (for recommended values see H&H pg. 72)
MR = 2 # Mixture Ratio

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
ablative = Ablative(engine, run_time=40, T_wall_i=3220, split=0)


# ======================= Steady-State Heat Transfer: Wall Thickness Sizing  ============================
t_A_design = 1.5*.0254 # ablator thickness (m) # if this if ablator was used as insulation
Max_OW_temp = 550 #K
liner_idx = 0 # uses the initial combustion chamber properties # This should be throat at some point
minT= 3200 # this is just to not have to do the entire temp regime
T_aw= ablative.size_ablator_thickness(burn_time=40,T_w=minT,t_A_design=t_A_design,wall_idx=liner_idx)[-3] # adiabatic wall temp is constant because we're looking at the throat only right now
maxT = T_aw
design_OW_thickness = .125*.0254 #m

# See overall trend by plotting a bunch of design wall temps:
Temp_arr_start = np.linspace(minT,maxT,100000)
t_OW_arr = []
T_OW_arr = []
Temp_arr=[]
Temp_arr2=[]
T_I_arr=[]
h_g_arr =[]

for i in Temp_arr_start:
    t_OW_i = ablative.size_ablator_thickness(burn_time=40,T_w=i,t_A_design=t_A_design,wall_idx=liner_idx)[1]
    T_I_i = ablative.size_ablator_thickness(burn_time=40,T_w=i,t_A_design=t_A_design,wall_idx=liner_idx)[-2]
    h_g_i = ablative.size_ablator_thickness(burn_time=40, T_w=i,t_A_design=t_A_design,wall_idx=liner_idx)[3]
    T_OW_i = ablative.size_ablator_thickness(burn_time=40, T_w=i,t_A_design=t_A_design,wall_idx=liner_idx)[-1]

    if T_I_i>0:
        Temp_arr2= np.append(Temp_arr2,i)
        T_I_arr = np.append(T_I_arr, T_I_i)

    if t_OW_i>0 and t_OW_i <(10*.0254):
        t_OW_arr =np.append(t_OW_arr,t_OW_i)
        T_OW_arr = np.append(T_OW_arr, T_OW_i)
        Temp_arr = np.append(Temp_arr,i)
        h_g_arr = np.append(h_g_arr,h_g_i)


'''plt.figure(1)
plt.scatter(Temp_arr,t_OW_arr/.0254,marker=".",label="wall Thickness (in)")
plt.scatter(T_aw*np.ones((len(Temp_arr),1)),t_OW_arr/.0254,label="Adiabatic Wall temp (K)",marker=".")
plt.scatter(Temp_arr,design_OW_thickness/.0254*np.ones((len(Temp_arr),1)),label="Design Wall Thickness (in)",marker=".")
plt.legend()
plt.xlabel("Design Wall Temperature (K)")
plt.ylabel("Calculated OW Wall Thickness")
plt.xlim([minT,maxT])
plt.show()


plt.figure(2)
plt.scatter(Temp_arr,T_OW_arr,label= "OW Temperature (K)",marker=".")
#plt.plot(Temp_arr,T_I_arr,label= "Interfacial Temperature (K)")
#plt.plot(T_aw*np.ones((len(Temp_arr),1)),T_I_arr,label="Adiabatic Wall temp")
plt.legend()
plt.xlabel("Design Wall Temperature (K)")
plt.ylabel("Temperature (K)")
plt.show()

plt.figure(3)
plt.scatter(Temp_arr2,T_I_arr,label= "Interfacial Temperature (K)",marker=".")
plt.legend()
plt.xlabel("Design Wall Temperature (K)")
plt.ylabel("Temperature (K)")
plt.show()'''

'''plt.figure(4)
plt.plot(Temp_arr,h_g_arr, label='Gas-Side Convective Coefficient')
plt.xlabel("Design Wall temp (K)")
plt.ylabel("Gas-Side Convective Coefficient (W/mK)")
plt.show()'''

# Pick individual points by selecting design wall temp (T_w_design):
# thickness here is the overwrap thickness
def size_thickness(ablative,OW_thickness):
    Temp_arr_start = np.linspace(minT, maxT, 1000)
    t_OW_arr = []
    T_OW_arr = []
    Temp_arr = []
    Temp_arr2 = []
    T_I_arr = []
    h_g_arr = []

    for i in Temp_arr_start:

        (_, t_OW_i, _, h_g_i, _, T_I_i, T_OW_i)= ablative.size_liner_thickness1(burn_time=40, T_w=i, wall_idx=liner_idx,
                                                                                t_A_design=t_A_design)

        if T_I_i > 0:
            Temp_arr2 = np.append(Temp_arr2, i)
            T_I_arr = np.append(T_I_arr, T_I_i)

        if t_OW_i > 0 and t_OW_i < (10 * .0254):
            t_OW_arr = np.append(t_OW_arr, t_OW_i)
            T_OW_arr = np.append(T_OW_arr, T_OW_i)
            Temp_arr = np.append(Temp_arr, i)
            h_g_arr = np.append(h_g_arr, h_g_i)


    idx = ablative.closest(t_OW_arr, OW_thickness)
    T_w = Temp_arr[idx]
    (t_A, t_OW, q_conv, h_g, T_aw, T_I, T_OW) = ablative.size_liner_thickness1(burn_time=40, T_w=T_w,
                                                                               wall_idx=liner_idx,
                                                                               t_A_design=t_A_design)
    return t_A, t_OW, q_conv, h_g, T_aw, T_I, T_OW,T_w


OW_thickness_arr =np.linspace(.05*.0254,.5*.0254,100)
T_I_arr = np.zeros((len(OW_thickness_arr),1))
T_OW_arr = np.zeros((len(OW_thickness_arr),1))
T_w_arr = np.zeros((len(OW_thickness_arr),1))
for i in range(len(OW_thickness_arr)):
    (t_A, t_OW, q_conv, h_g, T_aw, T_I, T_OW,T_w)= size_thickness(ablative,OW_thickness_arr[i])
    T_I_arr[i] = T_I
    T_OW_arr[i] = T_OW
    T_w_arr[i] = T_w

plt.figure()
plt.scatter(OW_thickness_arr/.0254,T_I_arr,label='Interface',marker=".")
plt.scatter(OW_thickness_arr/.0254,T_OW_arr,label='Outer Wall',marker=".")
plt.plot(OW_thickness_arr/.0254,np.ones((len(OW_thickness_arr),1))*Max_OW_temp,label="Max OW Temp")
plt.legend()
plt.xlabel("Outer Wall Design Thickness (in)")
plt.ylabel("Temperature (K)")
plt.show()

# CHOSEN DESIGN PARAMETERS
t_OW_CHOSEN = .125 *.0254 #m
DT_CHOSEN = 50 #K
(t_A, t_OW, q_conv, h_g, T_aw, T_I, T_OW,T_w)= size_thickness(ablative,t_OW_CHOSEN)

# OW Structural Analysis
mat_OW = ablative.CarbonEpoxy()
(stress_OW,FOS_OW) = ablative.chamber_stress_analysis(mat_OW,t_OW_CHOSEN,DT_CHOSEN)

print("")
print("===Chamber Sizing Analysis===")
print("Chosen Design OW Wall Thickness (in): "+str(t_OW_CHOSEN/.0254))
print("Design Wall Temperature(K): "+str(T_w))
print("Adiabatic Wall Temp AKA Hot-gas Temp at throat (K): " + str(T_aw))
print("Interface Temp (K): " + str(T_I))
print("Overwrap outer-wall temp: " + str(T_OW))
print("Heat flux rate: " + str(q_conv))
print("")
print("Ablator Thickness (in): " + str(t_A/.0254))
print("Overwrap Thickness (in): " + str(t_OW/.0254))
print("Gas Convective Coefficient (W/mK): " + str(h_g))
print("===Chamber Structural Analysis===")
print("Overwrap stress (MPa) " + str(stress_OW))
print("Overwrap FOS " + str(FOS_OW))
