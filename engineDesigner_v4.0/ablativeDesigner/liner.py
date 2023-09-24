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

# Daedelus

'''thrust = 19500 # Thrust [N]
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

ablative = Ablative(engine, run_time=40, T_wall_i=3220, split=0)


# ======================= Steady-State Heat Transfer: Wall Thickness Sizing  ============================
#t_A_design = 3*.0254 # ablator thickness (m)
liner_idx = 0 # uses the initial combustion chamber properties
t_A_design=[]
minT= 3200
T_aw= ablative.size_ablator_thickness(burn_time=40,T_w=minT,t_A_design=t_A_design,wall_idx=liner_idx)[-3] # adiabatic wall temp is constant because we're looking at the throat only right now
maxT = T_aw

design_wall_thickness = .125*.0254 #m

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





plt.figure(1)
plt.plot(Temp_arr,t_OW_arr/.0254,label="wall Thickness (in)")
plt.plot(T_aw*np.ones((len(Temp_arr),1)),t_OW_arr/.0254,label="Adiabatic Wall temp (K)")
plt.plot(Temp_arr,design_wall_thickness/.0254*np.ones((len(Temp_arr),1)),label="Design Wall Thickness (in)")
plt.legend()
plt.xlabel("Design Wall Temperature (K)")
plt.ylabel("Calculated OW Wall Thickness")
plt.xlim([minT,maxT])
plt.show()
idx = ablative.closest(t_OW_arr, design_wall_thickness)
T_w = Temp_arr[idx]




plt.figure(2)
plt.plot(Temp_arr,T_OW_arr,label= "OW Temperature (K)")
#plt.plot(Temp_arr,T_I_arr,label= "Interfacial Temperature (K)")
#plt.plot(T_aw*np.ones((len(Temp_arr),1)),T_I_arr,label="Adiabatic Wall temp")
plt.legend()
plt.xlabel("Design Wall Temperature (K)")
plt.ylabel("Temperature (K)")
plt.show()

plt.figure(3)
plt.plot(Temp_arr2,T_I_arr,label= "Interfacial Temperature (K)")
plt.legend()
plt.xlabel("Design Wall Temperature (K)")
plt.ylabel("Temperature (K)")
plt.show()

'''plt.figure(3)
plt.plot(Temp_arr,h_g_arr, label='Gas-Side Convective Coefficient')
plt.xlabel("Design Wall temp (K)")
plt.ylabel("Gas-Side Convective Coefficient (W/mK)")
plt.show()'''

# Pick individual points by selecting design wall temp (T_w_design):
# thickness here is the overwrap thickness
def size_thickness(ablative,thickness):
    idx = ablative.closest(t_OW_arr, thickness)
    T_w = Temp_arr[idx]
    (t_A, t_OW, q_conv, h_g, T_aw, T_I, T_OW) = ablative.size_ablator_thickness(burn_time=40, T_w=T_w,t_A_design=t_A_design,wall_idx=liner_idx)
    return t_A, t_OW, q_conv, h_g, T_aw, T_I, T_OW


thickness_arr =np.linspace(.05*.0254,.5*.0254,100)
T_I_arr = np.zeros((len(thickness_arr),1))
T_OW_arr = np.zeros((len(thickness_arr),1))

for i in range(len(thickness_arr)):
    (t_A, t_OW, q_conv, h_g, T_aw, T_I, T_OW)= size_thickness(ablative,thickness_arr[i])
    T_I_arr[i] = T_I
    T_OW_arr[i] = T_OW

plt.figure(4)
plt.plot(thickness_arr/.0254,T_I_arr,label='Interface')
plt.plot(thickness_arr/.0254,T_OW_arr,label='Outer Wall')
plt.legend()
plt.xlabel("Outer Wall Design Thickness (in)")
plt.ylabel(" Temperature (K)")
plt.show()


T_w_design = T_w
(t_A,t_OW,q_conv,h_g,T_aw,T_I,T_OW)=ablative.size_ablator_thickness(burn_time=40,T_w=T_w_design,t_A_design=t_A_design,wall_idx=liner_idx)
print("")
print("===Chamber Sizing Analysis===")
print("Design Wall Thickness (in): "+str(design_wall_thickness/.0254))
print("Design Wall Temperature(K): "+str(T_w))
print("Adiabatic Wall Temp AKA Hot-gas Temp at throat (K): " + str(T_aw))
print("Interface Temp (K): " + str(T_I))
print("Overwrap outer-wall temp: " + str(T_OW))
print("Heat flux rate: " + str(q_conv))
print("")
print("Ablator Thickness (in): " + str(t_A/.0254))
print("Overwrap Thickness (in): " + str(t_OW/.0254))
print("Gas Convective Coefficient (W/mK): " + str(h_g))

