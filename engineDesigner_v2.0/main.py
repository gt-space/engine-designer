# MAIN.PY – Master file for running full engine analysis
# Effectively the same as main.py in contourDesigner, but with regen functionality

from contourDesigner.design import Engine
from regenDesigner.regen import RegenJacket
import matplotlib.pyplot as plt

import numpy as np
from bartz import bartz


## GENERAL INPUT PARAMETERS ##
# For basic engine configuration and adjustment
thrust = 780 * 4.44822 # Thrust [lbf to N]
P_c = 260 * 0.0689476 # Chamber pressure at injector face [psi to bar]
P_e = 8 / 14.508 # Desired exit presure for expansion ratio [bar] (if in doubt put ambient)
con_rat = 4.82 # Contraction ratio
L_star = 1.05 # Characteristic length [m] (for recommended values see H&H pg. 72)
MR = 1.8# Mixture ratio by weight (ox/fuel)

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

wall_thickness = 0.002 # Thickness of liner (distance from inner wall to channel bottom) [m]
min_fin_width = 0.0005 # Fin width at the throat [m]
min_channel_width = 0.001 # Minimum allowable channel width [m] (determined by manufacturing capability)
var_channel_width = True# Turn variable channel width on or off (e.g., using an end mill vs a slitting saw)
starting_index = 0 # Index to begin analysis for non-full regen designs. Where the fuel enters


## CALL ENGINE DESIGN SCRIPTS ##
engine = Engine(thrust, P_c, P_e, con_rat, L_star = L_star, MR = MR, adv_data = adv_data) # Create engine object
engine.design_engine() # Run engine design process (source code in contourDesigner)

#np.savetxt('nozzle_contour_781_260_8.csv', engine.engineContour, delimiter=",")

## DISPLAY RESULTS ##
print(" +++ ENGINE DESIGN RESULTS +++")
print("Exit Velocity [m/s]: " + str(engine.V_exit))
print("Mass Flow Rate [kg/s]: " + str(engine.mDot_tot))
print("Fuel Mass Flow Rate [kg/s]: " + str(engine.mDot_f))
print("LOX Mass Flow Rate [kg/s]: " + str(engine.mDot_o))
print("Throat Area [m^2]: " + str(engine.A_t))
print("C-star [m/s]: " + str(engine.C_star))
print("Exit angle [deg]: " + str(engine.theta_e))
print("Chamber Radius [m]: " + str(engine.engineProps[0, 0]))
print("Exit Radius [m]: " + str(engine.engineProps[-1, 0]))
print(" ")

# Plot a specific value (change the 0 to the index of the property you want to see):
# Indecies can be found in design.py header
plt.figure()
plt.plot(engine.engineProps[:, 1], engine.engineProps[:, 0])
plt.xlabel('Distance from Injector [m]', fontsize=16)
plt.ylabel('Radius [m]', fontsize=16)
plt.xlim(0, .55)
plt.ylim(0, .55)
#plt.show()

#Optimize regen design
jacket = RegenJacket(engine, wall_thickness, min_fin_width, min_channel_width = min_channel_width, var_channel_width = var_channel_width, T_ci=300, start_ind = starting_index)# Create jacket object - regen_3
jacket.simulate_regen()

num_channels = 75
channel_w = .001

jacket.single_design(num_channels, channel_w)

## DISPLAY RESULTS ##
print(" +++ REGEN SIM RESULTS +++")
print("Number of Channels: " + str(num_channels))
print("Max wall temp [K]: " + str(max(jacket.wall_temps)[0]))
print("Pressure Loss [psi]: " + str((jacket.pressures[-1]-jacket.pressures[0])[0]*14.5038))
print("Coolant Outlet Temperature [K]: " + str(jacket.coolant_temps[0][0]))





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












