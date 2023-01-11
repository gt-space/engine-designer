#MAIN.PY – Master file for running full engine analysis
# Effectively the same as main.py in contourDesigner, but with regen functionality

from contourDesigner.design import Engine
from regenDesigner.regen import RegenJacket
import matplotlib.pyplot as plt

import numpy as np
from bartz import bartz
from chan_contour import get_chan_contour


## GENERAL INPUT PARAMETERS ##
# For basic engine configuration and adjustment
thrust = 3500 #781 * 4.44822 # Thrust [lbf to N]
P_c = 18 #261 * 0.0689476 # Chamber pressure at injector face [psi to bar]
P_e = 8 / 14.5038 # Desired exit presure for expansion ratio [bar] (if in doubt put ambient)
con_rat = 5.5 # Contraction ratio
L_star = 1.05 # Characteristic length [m] (for recommended values see H&H pg. 72)
MR = 1.8# Mixture ratio by weight (ox/fuel)
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
max_wall_temp = 1200

outlet_temp = 550
uts_FOS_min = 1.5


## CALL ENGINE DESIGN SCRIPTS ##
engine = Engine(thrust, P_c, P_e, con_rat, L_star = L_star, MR = MR, adv_data = adv_data, cstar_eff = cstar_eff) # Create engine object
engine.design_engine() # Run engine design process (source code in contourDesigner)

#Optimize regen design
jacket = RegenJacket(engine, wall_thickness, min_fin_width, min_channel_width = min_channel_width,
                     var_channel_width = var_channel_width,T_ci=300, start_ind = starting_index, max_wall_temp=max_wall_temp)# Create jacket object
MR = jacket.MR_optimization(outlet_temp, uts_FOS_min)

num_channels = 45
channel_w = 0.002
jacket.single_design(MR, num_channels, channel_w, plot=True)


##TOLERANCE ANALYSIS
#Channel height variation

d_c_h = 0.002
string = "c_h + " + str(d_c_h) + "\""
d_c_h = d_c_h * .0254
jacket1 = RegenJacket(engine, wall_thickness, min_fin_width, min_channel_width = min_channel_width,
                 var_channel_width = var_channel_width,T_ci=300, start_ind = starting_index, max_wall_temp=max_wall_temp)
jacket1.final_pass(MR, jacket.num_channels, jacket.channel_h_arr, jacket.channel_w_arr, jacket.wall_t)


jacket2 = RegenJacket(engine, wall_thickness, min_fin_width, min_channel_width = min_channel_width,
                 var_channel_width = var_channel_width,T_ci=300, start_ind = starting_index, max_wall_temp=max_wall_temp)
c_h_arr = jacket.channel_h_arr + d_c_h
jacket2.final_pass(MR, jacket.num_channels, c_h_arr, jacket.channel_w_arr, jacket.wall_t)
#T_wg_diff = jacket2.T_wg - jacket1.T_wg
#ult_FOS_diff = jacket2.ult_FOS_arr - jacket1.ult_FOS_arr


plt.figure()
plt.title("Inner Wall Temperature")
plt.xlabel("Distance from injector (m)")
plt.ylabel("Inner Wall Temp (K)")
plt.plot(jacket1.ax_pos, jacket1.wall_temps)
plt.plot(jacket2.ax_pos, jacket2.wall_temps)
plt.legend(["Nominal", string])

plt.figure()
plt.title("Coolant Temperature")
plt.xlabel("Distance from injector (m)")
plt.ylabel("Coolant Temp (K)")
plt.plot(jacket1.ax_pos, jacket1.coolant_temps)
plt.plot(jacket1.ax_pos, jacket2.coolant_temps)
plt.legend(["Nominal", string])

plt.figure()
plt.title("Ultimate FOS")
plt.xlabel("Distance from injector (m)")
plt.ylabel("Ultimate FOS")
plt.plot(jacket1.ax_pos, jacket1.ult_FOS_arr)
plt.plot(jacket2.ax_pos, jacket2.ult_FOS_arr)
plt.legend(["Nominal", string])

plt.show()






















## DISPLAY RESULTS ##
print(" +++ REGEN SIM RESULTS +++")
print("Number of Channels: " + str(num_channels))
print("Max wall temp [K]: " + str(max(jacket.wall_temps)[0]))
print("Pressure Loss [psi]: " + str((jacket.pressures[-1]-jacket.pressures[0])[0]*14.5038))
print("Coolant Outlet Temperature [K]: " + str(jacket.coolant_temps[0][0]))


def con_rat_fixer(num_channels, channel_w):
    con_rat_upper = 7
    con_rat_lower = 3

    con_rat_avg = (con_rat_upper + con_rat_lower) /2
    engine = Engine(thrust, P_c, P_e, con_rat_avg, L_star=L_star, MR=MR, adv_data=adv_data)  # Create engine object
    engine.design_engine()  # Run engine design process (source code in contourDesigner)
    jacket = RegenJacket(engine, wall_thickness, min_fin_width, min_channel_width=min_channel_width,
                         var_channel_width=var_channel_width, T_ci=300, start_ind=starting_index)
    jacket.single_design(num_channels, channel_w)
    rad_diff_avg = jacket.engine.engineContour[0, 0] + jacket.wall_t_arr[0] + jacket.channel_h_arr[0, 0] - (jacket.engine.engineContour[-1, 0] + jacket.wall_t_arr[-1] + jacket.channel_h_arr[-1, 0])


    while np.abs(rad_diff_avg) > .00000001:
        if rad_diff_avg > 0:
            con_rat_upper = con_rat_avg
        else:
            con_rat_lower = con_rat_avg

        con_rat_avg = (con_rat_upper + con_rat_lower) / 2

        engine = Engine(thrust, P_c, P_e, con_rat_avg, L_star=L_star, MR=MR, adv_data=adv_data)  # Create engine object
        engine.design_engine()  # Run engine design process (source code in contourDesigner)
        jacket = RegenJacket(engine, wall_thickness, min_fin_width, min_channel_width=min_channel_width,
                             var_channel_width=var_channel_width, T_ci=300, start_ind=starting_index)
        jacket.single_design(num_channels, channel_w)
        rad_diff_avg = jacket.engine.engineContour[0, 0] + jacket.wall_t_arr[0] + jacket.channel_h_arr[0, 0] - (jacket.engine.engineContour[-1, 0] + jacket.wall_t_arr[-1] + jacket.channel_h_arr[-1, 0])
    return con_rat_avg
#con_rat = con_rat_fixer(num_channels, channel_w)
#print("Contraction Ratio: " + str(con_rat))
'''
con_rat = 5.5
engine = Engine(thrust, P_c, P_e, con_rat, L_star=L_star, MR=MR, adv_data=adv_data)  # Create engine object
engine.design_engine()  # Run engine design process (source code in contourDesigner)
jacket = RegenJacket(engine, wall_thickness, min_fin_width, min_channel_width=min_channel_width,
                     var_channel_width=var_channel_width, T_ci=300, start_ind=starting_index)
jacket.single_design(num_channels, channel_w, plot=True)
get_chan_contour(jacket)
'''
def Isp_calc(): # ISP Calculation
    P_c_arr = np.arange(200, 801, 100) * .0689476
    MR_arr = np.arange(1.6, 2.6, .1)
    Isp_arr = np.empty([len(P_c_arr), len(MR_arr)])
    exit_area_arr = np.empty([len(P_c_arr), len(MR_arr)])
    for i in np.arange(len(P_c_arr)):
        for j in np.arange(len(MR_arr)):
            engine = Engine(thrust, P_c_arr[i], P_e, con_rat, L_star=L_star, MR=MR_arr[j], adv_data=adv_data)  # Create engine object
            engine.design_engine()  # Run engine design process (source code in contourDesigner)
            Isp_arr[i, j] = thrust / engine.mDot_tot / 9.8014
            exit_area_arr[i, j] = np.pi * engine.engineProps[-1, 0] ** 2

    out_arr = np.concatenate([Isp_arr, exit_area_arr], 0)
    np.savetxt('Isp_table.csv', out_arr, delimiter= ",")
    return

def thrust_MR_study():
    thrust_arr = np.arange(4000, 12001, 1000) * 4.44822
    MR_init_arr = [1.8, 2, 2.3, 2.5]
    MR_out_arr = np.empty([len(thrust_arr), len(MR_init_arr)])

    for i in np.arange(len(MR_init_arr)):
        engine = Engine(thrust_arr[0], P_c, P_e, con_rat, L_star=L_star, MR=MR_init_arr[i], adv_data=adv_data)
        engine.design_engine()
        station = 98  # np.argmin(engine.engineProps[:, 0])
        length = engine.engineProps[199, 1] - engine.engineProps[0, 1]
        heat_flux_allow = bartz(engine, 1000, station)[1] * (
                    2 * np.pi * engine.engineProps[station, 0]) * length / engine.mDot_f
        for j in np.arange(len(thrust_arr)):

            MR_upper = 2.6
            engine = Engine(thrust_arr[j], P_c, P_e, con_rat, L_star=L_star, MR=MR_upper, adv_data=adv_data)
            engine.design_engine()
            length = engine.engineProps[199, 1] - engine.engineProps[0, 1]
            heat_flux_diff_upper = bartz(engine, 1000, station)[1] * (
                        2 * np.pi * engine.engineProps[station, 0]) * length / engine.mDot_f - heat_flux_allow

            MR_lower = 1.4
            engine = Engine(thrust_arr[j], P_c, P_e, con_rat, L_star=L_star, MR=MR_lower, adv_data=adv_data)
            engine.design_engine()
            length = engine.engineProps[199, 1] - engine.engineProps[0, 1]
            heat_flux_diff_lower = bartz(engine, 1000, station)[1] * (
                        2 * np.pi * engine.engineProps[station, 0]) * length / engine.mDot_f - heat_flux_allow

            if heat_flux_diff_upper < 0:
                MR_avg = MR_upper
            elif heat_flux_diff_lower > 0:
                MR_avg = MR_lower
            else:

                MR_avg = (MR_upper + MR_lower) / 2
                engine = Engine(thrust_arr[j], P_c, P_e, con_rat, L_star=L_star, MR=MR_avg, adv_data=adv_data)
                engine.design_engine()
                length = engine.engineProps[199, 1] - engine.engineProps[0, 1]
                heat_flux_diff = bartz(engine, 1000, station)[1] * (
                            2 * np.pi * engine.engineProps[station, 0]) * length / engine.mDot_f - heat_flux_allow
                while np.abs(heat_flux_diff) > 10000:
                    if heat_flux_diff > 0:
                        MR_upper = MR_avg
                    else:
                        MR_lower = MR_avg
                    MR_avg = (MR_upper + MR_lower) / 2
                    engine = Engine(thrust_arr[j], P_c, P_e, con_rat, L_star=L_star, MR=MR_avg, adv_data=adv_data)
                    engine.design_engine()
                    length = engine.engineProps[199, 1] - engine.engineProps[0, 1]
                    heat_flux_diff = bartz(engine, 1000, station)[1] * (
                                2 * np.pi * engine.engineProps[station, 0]) * length / engine.mDot_f - heat_flux_allow
            MR_out_arr[j, i] = MR_avg
            print(length)

    plt.figure()
    plt.xlabel('Thrust (lbf)')
    plt.ylabel('MR')
    plt.title('Required MR holding q_in/mdot_f constant')
    plt.plot(thrust_arr / 4.44822, MR_out_arr[:, 0])
    plt.plot(thrust_arr / 4.44822, MR_out_arr[:, 1])
    plt.plot(thrust_arr / 4.44822, MR_out_arr[:, 2])
    plt.plot(thrust_arr / 4.44822, MR_out_arr[:, 3])
    plt.show()
    return
