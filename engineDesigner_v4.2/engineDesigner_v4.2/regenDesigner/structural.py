import numpy as np
import matplotlib.pyplot as plt
from .bartz import bartz


class CopperWall:

    def __init__(self):
        self.cte = 17.2 * 10 ** -6
        self.nu = .33


    def get_wall_props(self, T_wg):
        if T_wg > 1250:
            yts = 0
        elif T_wg > 805:
            yts = -31136377/(1250-805) * (T_wg - 805) + 31136377
        else:
            yts = -9.315 * 10 ** (-2) * T_wg ** 3 + 98.52 * T_wg ** 2 - 6.2521 * 10 ** 4 * T_wg + 66215000

        if T_wg > 1250:
            uts = 0
        elif T_wg > 1088.706:
            uts = -13961880/(1250-1088.706) * (T_wg - 1088.706) + 13961880
        elif T_wg > 805:
            uts = (13961880 - 57527820)/(1088.706-805) * (T_wg - 805) + 57527820
        else:
            uts = (-0.000001028 * (T_wg ** 3) + 0.001625753 * (T_wg ** 2) - 1.106140710 * (
            T_wg) + 430.709111917) * 10 ** 6



        return(uts, yts)

    def structural_analysis(self, delT, T_wg):
        axial_constraint = False
        radial_constraint = False

        # temp = jacket.wall_t_arr[i]

        (uts, yts) = self.get_wall_props(T_wg)

        E = (-0.0205 * T_wg + 121.22) * 10 ** 9


        K_533 = 90.9777 * 10 ** 6
        K_755 = 76.6298 * 10 ** 6
        K = K_533 * ((755.4 - T_wg) / (755.4 - 533.1)) + K_755 * (
                (T_wg - 533.1) / (755.4 - 533.1))  # Strength coeff, linear interp
        n_533 = 0.103928
        n_755 = 0.1261

        n = n_533 ** ((755.4 - T_wg) / (755.4 - 533.1)) * n_755 ** (
                (T_wg - 533.1) / (755.4 - 533.1))  # Strain hardening exp., log interp
        n = n_533 * ((755.4 - T_wg) / (755.4 - 533.1)) + n_755 * (
                (T_wg - 533.1) / (755.4 - 533.1))  # Strength coeff, linear interp
        delT_uniform = T_wg - 297
        # delT = T_wg - self.T_wc_arr[147]
        # Thermal expansion, uniform, unconstrained

        eps_ax = eps_rad = eps_circ = self.cte * delT_uniform

        eps_inner_equiv = 2 * self.cte * delT / (1 - self.nu)

        sig_inner_equiv = eps_inner_equiv * E

        if sig_inner_equiv > yts:
            sig_inner_equiv = K * eps_inner_equiv ** n

        if axial_constraint and radial_constraint:
            eps_inner_equiv = eps_inner_equiv * 2
            sig_inner_equiv = sig_inner_equiv * 2

        elif axial_constraint or radial_constraint:
            eps_inner_equiv = eps_inner_equiv * 1.6
            sig_inner_equiv = sig_inner_equiv * 1.6

        eps_min = self.cte * (90 - 297)
        tot_strain_range = eps_inner_equiv - eps_min

        FOS_yield = yts / sig_inner_equiv
        FOS_ult = uts / sig_inner_equiv

        # fin_w = ((2 * np.pi * (jacket.engine.engineContour[i, 0] + jacket.wall_t_arr[i])) - (jacket.num_channels * jacket.channel_w_arr[i])) / jacket.num_channels
        # m = np.sqrt(jacket.h_c_arr[i] * 2 / (jacket.cond_w * fin_w))

        # thet_rat = 1 / np.cosh(m * jacket.channel_h_arr[i])
        # temp_fintip = thet_rat * (jacket.T_wc_arr[i] - jacket.coolant_temps[i]) + jacket.coolant_temps[i]
        return (sig_inner_equiv, FOS_yield, FOS_ult)
copper = CopperWall()
(sig_inner_equiv, FOS_yield, FOS_ult) = copper.structural_analysis(13, 823.6)
print(sig_inner_equiv)
print(FOS_yield)
print(FOS_ult)
'''
    def struct_optimization(self, i, plot=False):
        A_g = np.pi * self.engine.engineProps[i, 0] ** 2
        T_wg_arr = np.arange(300, 925, 10)
        FOS_y_arr = np.zeros(len(T_wg_arr))
        FOS_ult_arr = np.zeros(len(T_wg_arr))
        sig_arr = np.zeros(len(T_wg_arr))
        uts_arr = np.zeros(len(T_wg_arr))
        yts_arr = np.zeros(len(T_wg_arr))
        for j in np.arange(len(T_wg_arr)):
            (h_g, q_in, T_aw) = bartz(self.engine, T_wg_arr[j], i)
            delT = q_in / (self.cond_w / self.wall_t_arr[i])
    
            (sig_arr[j], FOS_y_arr[j], FOS_ult_arr[j], uts_arr[j], yts_arr[j]) = self.structural_analysis(delT, T_wg_arr[j])
        if plot:
            plt.figure()
            plt.plot(T_wg_arr, FOS_ult_arr)
            plt.xlabel('Inner Wall Temp, K')
            plt.ylabel('FOS, Ultimate')
            plt.title('Ultimate Factor of Safety')
    
            plt.figure()
            plt.plot(T_wg_arr, FOS_y_arr)
            plt.xlabel('Inner Wall Temp, K')
            plt.ylabel('FOS, Yield')
            plt.title('Yield Factor of Safety')
    
            plt.figure()
            plt.plot(T_wg_arr, sig_arr)
            plt.plot(T_wg_arr, uts_arr)
            plt.plot(T_wg_arr, yts_arr)
            plt.xlabel('Inner Wall Temp, K')
            plt.ylabel('Stress, MPA')
            plt.title('Inner Wall Equivalent Stress')
            plt.legend(["Engine Stress", "Ultimate Stress", "Yield Stress"])
            plt.show()
        return
'''
