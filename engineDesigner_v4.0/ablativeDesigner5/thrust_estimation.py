import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve
from design import Engine
from bartz import bartz


def thrust_curve(self, throat_insert_ablation_rate, chamber_ablation_rate, burnTime, debugMach, showplot=True):
    burnTimes = np.arange(0, burnTime, .01)
    Thrust = np.zeros(burnTimes.size)
    Thrust_mdot_cnst = np.zeros(burnTimes.size)
    M_e = np.zeros(burnTimes.size)
    for i in range(len(burnTimes)):
        (Thrust_mdot_cnst[i], Thrust[i], M_e[i]) = self.thrust_change_throat(throat_insert_ablation_rate,
                                                                             chamber_ablation_rate, burnTimes[i])
    # Output the data

    '''
    print(" +++ THRUST GAIN CURVE ESTIMATION for LINER material +++")
    ablative = Ablative(engine, run_time=10,T_wall_i=0, split=0)
    time = 40 #seconds
    throat_insert_ablation_rate = 0.000508 # m/s # From Sutton table 14-5
    new_thrust = ablative.thrust_change_throat(throat_insert_ablation_rate, time)[0]
    ablative.thrust_curve(throat_insert_ablation_rate, time,debugMach=False,showplot=False)

    print("OLD THRUST: " + str(thrust))
    print("NEW FINAL THRUST: " + str(new_thrust))
    print("Thrust change: " + str(new_thrust-thrust) + " (N) " +str((new_thrust-thrust)*0.224809) + " (lbf), " +str((new_thrust-thrust)/thrust*100) + " %")

    plt.savefig('ThrustGainLINER.png')

    print(" +++ THRUST GAIN CURVE ESTIMATION for THROAT INSERT material +++")
    ablative = Ablative(engine, run_time=10,T_wall_i=0, split=0)
    time = 40 #seconds
    throat_insert_ablation_rate = 0.0001524 # m/s # From Sutton table 14-5
    new_thrust = ablative.thrust_change_throat(throat_insert_ablation_rate, time)[0]
    ablative.thrust_loss_curve(throat_insert_ablation_rate, time,debugMach=False,showplot=False)
    plt.savefig('ThrustGainThroatInsert.png')
    print("OLD THRUST: " + str(thrust))
    print("NEW FINAL THRUST: " + str(new_thrust))
    print("Thrust change: " + str(new_thrust-thrust) + " (N) " +str((new_thrust-thrust)*0.224809) + " (lbf), " +str((new_thrust-thrust)/thrust*100) + " %")


    '''
    # Plot the data
    if showplot:
        plt.figure()
        plt.plot(burnTimes, Thrust, 'b*', label='Thrust with changing Mdot')
        plt.plot(burnTimes, Thrust_mdot_cnst, 'g*', label='Thrust w/o changing Mdot')
        # Set the title and axis labels
        plt.title("Ablative Engine Thrust")
        plt.xlabel("Burn Time (s)")
        plt.ylabel("Thrust")
        plt.legend()
        # Show the plot
        plt.show()

    if debugMach:
        # Plot the data
        plt.figure()
        plt.plot(burnTimes, M_arr, 'b*')

        # Set the title and axis labels
        plt.title("Ablative Engine Exit Mach")
        plt.xlabel("Burn Time (s)")
        plt.ylabel("M_e")

        # Show the plot
        plt.show()
    return


def thrust_change_throat(self, throat_insert_ablation_rate, chamber_ablation_rate, time):
    Abl_rate_thr = throat_insert_ablation_rate
    # Engine properties
    A_e = np.pi * self.r_e ** 2

    # Eroded throat properties
    r_t = self.D_t / 2
    eroded_r_t = r_t - Abl_rate_thr * time
    eroded_A_t = np.pi * eroded_r_t ** 2
    eroded_epsilon = A_e / eroded_A_t

    # Eroded engine properties

    M_e = self.AonAStar_mach(eroded_epsilon, self.gamma_e)[1]

    p_0_over_p, T_0_over_T, rho_0_over_rho, A_over_Astar, MFP = self.compflowtool(self.gamma_e, M_e)
    T_e = T_0_over_T ** -1 * self.T_0
    R_e = 8314.5 / self.MW_e
    u_e = M_e * np.sqrt(self.gamma_e * R_e * T_e)

    # if mass flow rate doesnt change:
    Thrust_mdot_cnst = self.mdot * u_e

    # if mass flow rate DOES change as a result of chamber pressure dropping:

    gamma_c = self.engine.engineProps[0, 15]
    gamma = (gamma_c + self.gamma_e) / 2
    MW_c = self.engine.engineProps[0, 13]
    MW = (MW_c + self.MW_e) / 2
    R = 8314.5 / (MW)
    p_c_new = self.p_c_drop(chamber_ablation_rate, time, showplot=False)
    A = p_c_new / np.sqrt(R * self.T_0) * eroded_A_t
    B = np.sqrt(gamma * (2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1)))
    newmdot = A * B
    Thrust = newmdot * u_e
    # This assumes the following:
    # Chamber temperature does not change due to ablation
    # Molecular Weight (MW) and specific heat ratio (gamma) are average of the combustion chamber value and exit value
    return Thrust_mdot_cnst, Thrust, M_e


def p_c_drop(self, chamber_ablation_rate, burnTime, showplot=True):
    erate = chamber_ablation_rate  # m/s
    burnTimes = np.arange(0, burnTime, .01)  # SECONDS
    p1 = self.engine.P_inj * 1e5  # Pascals!
    ri = self.engine.engineProps[0, 0]  # Meters!
    l = self.engine.engineProps[self.engine.throatInd, 1]  # Meters

    V1 = np.pi * ri ** 2 * l
    r_final = ri + (erate * burnTime)
    Vf = np.pi * (r_final) ** 2 * l
    pf = p1 * V1 / Vf

    if showplot:
        p2 = []
        burnTimes = np.arange(0, burnTime, .01)
        for i in burnTimes:
            rnew = ri + (erate * i)
            V2 = np.pi * (rnew) ** 2 * l
            p2_val = p1 * V1 / V2
            p2 = np.append(p2, p2_val)

        t = burnTimes
        name = ['Erate: ' + str(erate)]
        plt.figure()
        plt.show()
        plt.plot(t, p2 / 6894.76, label=name)
        plt.xlabel('Time (s)')
        plt.ylabel('Chamber Pressure (psi)')

    return pf
