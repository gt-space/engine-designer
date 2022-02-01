import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from fuel_props import get_props
from bartz import bartz
from bartz import bartz_sigma
from scipy.optimize import fsolve
from scipy.signal import savgol_filter
from tqdm import tqdm

import sympy as sym

np.set_printoptions(threshold=sys.maxsize)  # Print full arrays (for debugging)


class RegenJacket:
    # Initialize the jacket object upon declaration
    def __init__(self, engine, wall_t=0.002, min_fin_w=0.001, min_channel_width=.001, var_channel_width=False, T_ci=300,
                 start_ind=0):
        self.engine = engine  # Engine object to be jacketed

        self.wall_t = wall_t  # Inner wall thickness [m]
        self.min_fin_w = min_fin_w  # Minimum desired fin width [m]
        self.min_channel_w = min_channel_width
        self.var_channel_w = var_channel_width
        self.T_ci = T_ci  # Coolant inlet temperature [K]
        self.start_ind = start_ind  # Index where fuel enters (0 is exit plane)

        self.numPTS = np.size(self.engine.engineProps, 0)  # Number of engine stations
        max_wall_temp = 1000  # Add in FEA/structural calcs to fix this later
        self.max_wall_temp_arr = np.ones([self.numPTS - 1 - self.start_ind, 1]) * max_wall_temp

        self.cond_w = 350

    def init_channel(self):
        #Initialize the channel, called each time a new possible design is started
        self.T_cb = self.T_ci  # Initial bulk temp [K]
        self.P_co = self.engine.P_inj * 1.2  # Initial coolant pressure [bar] 1.2 included as stiffness

        self.wall_temps = np.zeros([self.numPTS - 1 - self.start_ind, 1])  # Wall temperature vector
        self.coolant_temps = np.zeros([self.numPTS - 1 - self.start_ind, 1])  # Coolant temperature vector
        self.pressures = np.zeros([self.numPTS - 1 - self.start_ind, 1])  # Coolant pressure vector
        self.h_g_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1])  # Hot Gas side Convection coefficient vector
        self.h_c_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1])  # Coolant side Convection coefficient vector
        self.T_wc_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1])  # Coolant side wall temp vector
        self.channel_h_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1])  # Channel height array
        self.channel_w_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1])  # Channel width array
        self.fin_w_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1])  # Fin width array
        self.dP_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1])  # Pressure loss array
        self.Re_c_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1])  # Coolant Reynolds number array
        self.q_in_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1])  # Heat influx array
        self.q_out_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1]) # Heat outflux array
        self.A_effc_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1]) # Effective cooling area array

    # Friction correlation function
    def get_friction(self, Re, D_hyd):
        def colebrook(Re, D_hyd):

            k = 0.005 * (10 ** (-3))  # Surface roughness, guess based on our manufacuting capability, mm => m
            #Currently using Halbrook correlation to boost computational speed (generally accurate to ~3%)
            #To improve accuracy at significant computational expense, can solve numerically using bisection method
            '''
            fric_upper = 0.25
            fric_lower = 0.01
            fric_avg = 0.05
            fric_diff = -2*np.log10(k/(3.71*D_hyd) + 2.51/(Re*np.sqrt(fric_avg))) - 1.0/np.sqrt(fric_avg)
            i = 0
            while np.abs(fric_diff)  > .00001 and i < 50:
                if fric_diff > 0:
                    fric_upper = fric_avg
                else:
                    fric_lower = fric_avg
                fric_avg = (fric_upper + fric_lower) / 2
                fric_diff = -2 * np.log10(k / (3.71 * D_hyd) + 2.51 / (Re * np.sqrt(fric_avg))) - 1.0 / np.sqrt(fric_avg)
                i += 1


            fric = fric_avg
            if i == 50:
                fric = 1 / ((-1.8 * np.log10((k / D_hyd / 3.7) ** 1.11 + 6.9 / Re)) ** 2)
            #func = lambda x: -2*np.log10(k/(3.71*D_hyd) + 2.51/(Re*np.sqrt(x))) - 1.0/np.sqrt(x)
            #fric = fsolve(func, 1 / ((-1.8*np.log10((k / D_hyd / 3.7)**1.11 + 6.9/Re))**2))
            '''
            fric = 1 / ((-1.8 * np.log10((k / D_hyd / 3.7) ** 1.11 + 6.9 / Re)) ** 2)
            return fric

        # Solve for Darcy Friction Factor Approximation using Reynold's number
        if Re <= 2320:
            fric = 64 / Re  # Laminar flow approximation
        else:
            fric = colebrook(Re, D_hyd)  # Colebrook-White Equation solution
        return fric

    def qout_calc(self, channel_h, channel_w, i):
        '''
        Computes necessary parameters for coolant heat flux calculation, as well as local channel pressure loss
        Inputs: channel height, channel width, current station
        Outputs: coolant convection coefficient, effective cooling area, pressure loss
        '''
        (rho_c, C_pc, cond_c, viscK_c) = get_props(self.T_cb)
        R = self.engine.engineProps[i, 0]  # Radius [m]
        length = self.engine.engineProps[i, 1] - self.engine.engineProps[i - 1, 1]  # Station Length [m]
        length = math.sqrt(length ** 2 + (R - self.engine.engineProps[i - 1, 0]) ** 2)  # Account for angle at station
        fin_w = ((2 * math.pi * (R + self.wall_t)) - (self.num_channels * channel_w)) / self.num_channels  # fin width
        D_hyd = 2 * channel_w * channel_h / (channel_w + channel_h)  # Channel hydraulic diameter
        A_cc = channel_w * channel_h  # Cross sectional area of channel (m^2)
        vel_c = self.mDot_chan / (A_cc * rho_c)  # Coolant velocity in channel (m/s)
        Re_c = (vel_c * D_hyd) / viscK_c  # Coolant Reynolds number
        Pr_c = viscK_c * rho_c * C_pc / cond_c  # Coolant Prandtl Number

        # Apply Nusselt relation (see H&H pg 90):
        # h_c = 0.021 * (Re_c)**0.8 * (Pr_c)**0.4 * (0.64 + 0.36 * (T_cb/T_wc)) * (cond_c/D_hyd) # Old code
        # viscK_c_w = get_props(T_wc)[3]  # Get viscosity at wall
        # h_c = 0.0214 * Re_c**0.8 * (Pr_c)**0.4 * (1000000*viscK_c/(0.2))**0.14 * (cond_c/D_hyd)
        # h_c_2 = 0.021 * (Re_c**0.8) * (Pr_c**0.4) * (0.64+0.36*(T_cb/T_wc)) * (cond_c/D_hyd)
        # h_c = 0.023 * (Re_c ** 0.8) * (Pr_c ** 0.4) * (cond_c / D_hyd)  # Dittus Boelter reccomended for small dT Convection

        # Friction correlation

        f = self.get_friction(Re_c, D_hyd)
        Nu = (f / 8) * (Re_c - 1000) * Pr_c / (1 + 12.7 * ((f / 8) ** 0.5) * (Pr_c ** (2 / 3) - 1))
        h_c = Nu * cond_c / D_hyd * (1000000 * viscK_c / 0.2) ** 0.11

        # Future versions should edit this out of calculation
        # Run fin efficiency calculations to find the effective contact area
        # m = math.sqrt((h_c * perim_conf)/(cond_w * A_conf)) # Old fin efficiency parameter
        m = math.sqrt((2 * h_c) / (self.cond_w * fin_w))  # Fin efficiency length independent
        eta_fin = math.tanh(m * channel_h) / (m * channel_h)  # Fin efficiency
        A_fin = 2 * channel_h * length  # Area of single fin in analysis segment
        A_finTot = self.num_channels * A_fin  # Total fin area
        A_wallTot = self.num_channels * channel_w * length  # Total wall contact area
        A_totc = A_finTot + A_wallTot  # Total coolant side heat transfer area
        eta_tot = (A_wallTot + A_finTot * eta_fin) / A_totc  # Overall convection efficiency

        A_effc = A_totc * eta_tot

        dP = f * (length / D_hyd) * rho_c * (
                    (vel_c ** 2) / 2) / 100000  # Solve for final dP across passage segment (Pa > Bar)
        return h_c, A_effc, dP

    def height_optimization(self, channel_w, i, final=False):
        '''
        Computes the require channel height to close the thermal design for a given channel width
        Inputs: channel width, current station
        Outputs: station pressure loss, station coolant temp increase (channel height stored as jacket attribute
        '''

        # Initialize parameters
        (rho_c, C_pc, cond_c, viscK_c) = get_props(self.T_cb)  # Get coolant fluid properties

        R = self.engine.engineProps[i, 0]  # Radius [m]
        length = self.engine.engineProps[i, 1] - self.engine.engineProps[i - 1, 1]  # Station Length [m]
        length = math.sqrt(length ** 2 + (R - self.engine.engineProps[i - 1, 0]) ** 2)  # Account for angle at station

        T_wg = self.max_wall_temp_arr[i - 1]  # DESIRED wall temp
        (h_g, q_in, T_aw) = bartz(self.engine, T_wg, i)  # Gas side heat flux parameters
        A_totg = 2 * math.pi * R * length  # total hot gas side wall area
        q_in = q_in * A_totg / length  # Heat flux to heat transfer (length independent bc some of the lengths are funky)
        T_wc = T_wg - (q_in / (self.cond_w / self.wall_t) / (A_totg / length))  # Fix later to add in variable wall cooling

        fin_w = ((2 * math.pi * (R + self.wall_t)) - (self.num_channels * channel_w)) / self.num_channels  # fin width
        Re_c_min = 3000  # Minimum allowable Reynolds number in coolant, to avoid transition region

        channel_h_max_Re = 2 * self.mDot_chan / (Re_c_min * rho_c * viscK_c) - channel_w  # Channel height upper bound, use as initial value for optimization
        channel_h_max_buckling = (2 * np.pi * R - self.num_channels * channel_w) / self.num_channels * 3  # max fin height before buckling, as a multiple of fin width (FIX THIS LATER)
        channel_h_max_manuf = .75 * .0254  # max channel height due to radius of slitting saw

        channel_h_1 = min([channel_h_max_Re, channel_h_max_buckling, channel_h_max_manuf])  # Max channel height
        channel_h_2 = .0001  # Min channel height

        if channel_h_1 < channel_h_2:

            return (float('inf'), float('inf'))

        # q_in - q_out at max channel height
        (h_c, A_effc, dP) = self.qout_calc(channel_h_1, channel_w, i)
        q1 = h_c * A_effc / length * (T_wc - self.T_cb) - q_in

        # q_in - q_out at min channel height
        (h_c, A_effc, dP) = self.qout_calc(channel_h_2, channel_w, i)
        q2 = h_c * A_effc / length * (T_wc - self.T_cb) - q_in


        if q2 < 0: #q_out - q_in < 0 @ min channel height -> Invalid design

            return(float('inf'), float('inf'))

        elif q1 > 1000: #q_out - q_in > 0 @
            (h_c, A_effc, dP) = self.qout_calc(channel_h_1, channel_w, i)

            (tau, T_aw) = bartz_sigma(self.engine, i)

            T_0g = self.engine.engineProps[0, 9]  # Stagnation Temp
            k_t = self.cond_w / self.wall_t
            gam = self.engine.engineProps[i, 15]  # Index gamma; No conversion needed
            mach = self.engine.engineProps[i, 4]
            a = (1 + ((gam - 1) / 2) * mach ** 2)
            func = lambda T_wg: T_wg * (tau * ((.5 * T_wg / T_0g * a + .5) ** -.68) * (a ** -.12)) * k_t * (
                        A_totg ** 2) - \
                                T_aw * (tau * ((.5 * T_wg / T_0g * a + .5) ** -.68) * (a ** -.12)) * k_t * (
                                            A_totg ** 2) + \
                                T_wg * (tau * ((.5 * T_wg / T_0g * a + .5) ** -.68) * (a ** -.12)) * h_c * (
                                            A_totg * A_effc) + \
                                T_wg * h_c * k_t * (A_totg * A_effc) - \
                                T_aw * (tau * ((.5 * T_wg / T_0g * a + .5) ** -.68) * (a ** -.12)) * h_c * (
                                            A_totg * A_effc) - \
                                h_c * self.T_cb * k_t * (A_totg * A_effc)  # ouch
            T_wg = fsolve(func, 1000)
            (h_g, q_in, T_aw) = bartz(self.engine, T_wg, i)  # Gas side heat flux parameters
            q_in = h_g * A_totg / length * (T_aw - T_wg)
            T_wc = T_wg - (q_in * length / k_t / A_totg)
            q_out = h_c * A_effc / length * (T_wc - self.T_cb)

        else: #q_out = q_in somewhere within allowable channel range

            channel_h_avg = (channel_h_2 + channel_h_1) / 2
            (h_c, A_effc, dP) = self.qout_calc(channel_h_avg, channel_w, i)
            q_avg = h_c * A_effc / length * (T_wc - self.T_cb) - q_in
            while np.abs(q_avg) > 1000:
                if q_avg > 0:
                    channel_h_2 = channel_h_avg
                else:
                    channel_h_1 = channel_h_avg
                channel_h_avg = (channel_h_1 + channel_h_2) / 2
                (h_c, A_effc, dP) = self.qout_calc(channel_h_avg, channel_w, i)
                q_avg = h_c * A_effc / length * (T_wc - self.T_cb) - q_in

            q1 = q_avg
            channel_h_1 = channel_h_avg
            q_out = q1 + q_in
            '''
            channel_h_1 = channel_h_2
            (h_c, A_effc, dP) = self.qout_calc(channel_h_1, channel_w, i)
            q_avg = h_c * A_effc / length * (T_wc - self.T_cb) - q_in

            while np.abs(q_avg) > 1000:

                (h_c, A_effc, dP) = self.qout_calc(channel_h_1+.000001, channel_w, i)
                q_avg2 = h_c * A_effc / length * (T_wc - self.T_cb) - q_in
                deriv = (q_avg2 - q_avg) / .000001

                channel_h_1 = channel_h_1 - q_avg/deriv

                (h_c, A_effc, dP) = self.qout_calc(channel_h_1, channel_w, i)
                q_avg = h_c * A_effc / length * (T_wc - self.T_cb) - q_in


            q_out = q_avg + q_in
            '''
        channel_h = channel_h_1

        dT_c = (q_in + q_out) * length / 2 / (C_pc * self.engine.mDot_f)  # Temperature change in coolant

        if (not self.var_channel_w) or (self.var_channel_w and final):
            self.T_cb += dT_c
            self.P_co -= dP
            self.channel_h_arr[i - 1, 0] = channel_h
            self.channel_w_arr[i - 1, 0] = channel_w
            self.fin_w_arr[i - 1, 0] = fin_w
            self.wall_temps[i - 1, 0] = T_wg
            self.coolant_temps[i - 1, 0] = self.T_cb
            self.pressures[i - 1, 0] = self.P_co
            self.h_g_arr[i - 1, 0] = h_g
            self.h_c_arr[i - 1, 0] = h_c
            self.T_wc_arr[i - 1, 0] = T_wc
            self.dP_arr[i - 1, 0] = dP
            self.Re_c_arr[i - 1, 0] = self.mDot_chan / (rho_c * viscK_c) * 2 / (
                        channel_w + channel_h)  # Double check that this is right
            self.q_in_arr[i - 1, 0] = q_in
            self.q_out_arr[i - 1, 0] = q_out
            self.A_effc_arr[i - 1, 0] = A_effc / length

        return (dP, dT_c)

    def fixed_width_optimization(self):
        '''
        Primary script for optimizing channel geometry assuming FIXED width channels
        Assumes a slitting saw is used to cut channels so width variation is not allowed along contour

        Loops through possible numbers of channels based on throat circumference:
            Loops through all possible channel widths for given num channels
                Loops through all stations in engine (bottom -> top) given channel width & num channels
                    Computes channel height at each station
                    Computes dP and coolant temp increase for given design

        Displays 3D Plots of dP and coolant temp increase as a function of number of channels and channel width
        '''
        print('PERFORMING FIXED WIDTH, VARIABLE HEIGHT CHANNEL OPTIMIZATION')
        max_channels = 2 * np.pi * np.min(self.engine.engineProps[:, 0]) / (self.min_channel_w + self.min_fin_w)  # Maximum allowable number of channels, constrained by throat
        num_channels_arr = np.arange(40, max_channels,10)  # Array of possible number of channels to iterate through (lower bound changed manually)
        a = 2 * np.pi * np.min(self.engine.engineProps[:, 0]) / np.min(num_channels_arr)  # Filler variable for max channel width
        channel_w_range_max = np.arange(0, a + .0001, .0001)  # Array of overall possible channel widths to iterate through (e.g., for smallest num of channels)
        outlet_temp_3d = np.ones([len(num_channels_arr), len(channel_w_range_max)]) * 550  # Coolant outlet temp for all possibilities #Edit 550 to make graph look pretty
        dP_tot_3d = np.zeros([len(num_channels_arr), len(channel_w_range_max)])  # dP array for all possibilities
        pbar_num_channels = tqdm(desc="Iterating through Possible Numbers of Channels", total=len(num_channels_arr), leave=False)
        for j in np.arange(len(num_channels_arr)):
            self.num_channels = num_channels_arr[j]
            self.mDot_chan = self.engine.mDot_f / self.num_channels

            channel_w_high = ((np.min(self.engine.engineProps[:,0]) + self.wall_t) * np.pi * 2 - self.min_fin_w * self.num_channels) / self.num_channels  # Max channel width
            channel_w_low = self.min_channel_w  # Min channel width
            channel_w_range = np.arange(channel_w_low, channel_w_high, .0001)  # Channel width array
            dP_tot_arr = np.zeros(len(channel_w_range))  # dP array for each possible channel width
            outlet_temp_arr = np.zeros(len(channel_w_range))  # outlet temp array for each possible channel width
            pbar_width = tqdm(desc="Iterating through Possible Channel Widths", total=len(channel_w_range), leave=False)
            for k in np.arange(len(channel_w_range)):
                self.init_channel()  # Initialize channel - resets everything back to initial conditions from previous iterations

                x = 0
                pbar_station = tqdm(desc="Optimizing Height Along Nozzle Contour",
                                    total=self.numPTS - 1 - self.start_ind, leave=False)
                while x < (self.numPTS - 1 - self.start_ind):
                    i = self.numPTS - x - self.start_ind - 1  # Reverse order
                    self.height_optimization(channel_w_range[k], i)  # Optimize channel height for given width
                    pbar_station.update(1)
                    x += 1
                pbar_station.close()
                outlet_temp_arr[k] = self.coolant_temps[0]
                dP_tot_arr[k] = np.sum(self.dP_arr)

                pbar_width.update(1)
            pbar_width.close()

            # Insert outlet_temp_arr and dP_tot_arr into proper position in overall arrays
            idx_in = int(self.min_channel_w / .0001)
            idx_out = int(idx_in + len(outlet_temp_arr))
            outlet_temp_3d[j, idx_in:idx_out] = outlet_temp_arr
            dP_tot_3d[j, idx_in:idx_out] = dP_tot_arr

            pbar_num_channels.update(1)
        pbar_num_channels.close()

        (channel_w_mesh, num_channels_mesh) = np.meshgrid(channel_w_range_max, num_channels_arr)

        # Plot data
        plt.figure()
        ax_ch = plt.axes(projection='3d')
        ax_ch.plot_surface(channel_w_mesh, num_channels_mesh, outlet_temp_3d, rstride=1, cstride=1,
                           cmap='viridis', edgecolor='none')
        ax_ch.set_xlabel('Channel Width (m)')
        ax_ch.set_ylabel('Number of Channels')
        ax_ch.set_zlabel('Temp (K)')
        plt.title('Coolant Outlet Temperature')

        plt.figure()
        ax_ch = plt.axes(projection='3d')
        ax_ch.plot_surface(channel_w_mesh, num_channels_mesh, dP_tot_3d, rstride=1, cstride=1,
                           cmap='viridis', edgecolor='none')
        ax_ch.set_xlabel('Channel Width (m)')
        ax_ch.set_ylabel('Number of Channels')
        ax_ch.set_zlabel('dP (bar)')
        plt.title('Channel dP')
        plt.show()

        return (outlet_temp_3d, dP_tot_3d)

    def var_width_optimization(self):
        '''
        Primary script for optimizing channel geometry assuming VARIABLE width channels
        Assumes an endmill is used to cut channels so width variation is allowed along contour

        Loops through possible numbers of channels based on throat circumference:
            Loops through all stations in engine (bottom -> top) given num channels
                Computes optimal channel height and width at each station by choosing the combination which minimizes dP while staying within 2% of minimum dT_c
                    (^This methodology should be improved)
                    Computes dP and coolant temp increase for given design

        Displays 2D Plots of dP and coolant temp increase as a function of number of channels
        '''
        print('PERFORMING VARIABLE WIDTH AND HEIGHT CHANNEL OPTIMIZATION')
        max_channels = 2 * np.pi * np.min(self.engine.engineProps[:, 0]) / (self.min_channel_w + self.min_fin_w)
        num_channels_arr = np.arange(50, max_channels, 10)
        dP_arr_tot = np.zeros(len(num_channels_arr))
        outlet_temp_arr_tot = np.zeros(len(num_channels_arr))
        pbar_num_channels = tqdm(desc="Iterating through Possible Numbers of Channels", total=len(num_channels_arr),
                                 leave=False)
        for j in np.arange(len(num_channels_arr)):
            self.num_channels = num_channels_arr[j]
            self.mDot_chan = self.engine.mDot_f / self.num_channels

            self.init_channel()

            pbar_station = tqdm(desc="Optimizing Channel Width/Height along Nozzle Contour",
                                total=self.numPTS - 1 - self.start_ind, leave=False)
            x = 0
            while x < (self.numPTS - 1 - self.start_ind):
                i = self.numPTS - x - self.start_ind - 1  # Reverse order
                channel_w_high = ((self.engine.engineProps[i, 0] + self.wall_t) * np.pi * 2 - self.min_fin_w * self.num_channels) / self.num_channels #Max channel width
                channel_w_low = self.min_channel_w #Min Channel width
                channel_w_range = np.arange(channel_w_low, channel_w_high, .0001)

                dP_arr = np.zeros(len(channel_w_range))
                dT_c_arr = np.zeros(len(channel_w_range))

                #Sort through coolant temp increases, only keep those which are within 2% of min value
                for k in np.arange(len(channel_w_range)):
                    (dP_arr[k], dT_c_arr[k]) = self.height_optimization(channel_w_range[k], i)
                a = np.min(dT_c_arr) * 1.02
                for l in np.arange(len(dT_c_arr)):
                    if dT_c_arr[l] > a:
                        dP_arr[l] = float('inf')
                #Choose min value in dP array (which is also within 2% of min coolant temp increase)
                a = np.argmin(dP_arr)

                self.height_optimization(channel_w_range[a], i, final=True)
                pbar_station.update(1)
                x += 1
            dP_arr_tot[j] = np.sum(self.dP_arr)
            outlet_temp_arr_tot[j] = self.coolant_temps[0]
            pbar_station.close()
            pbar_num_channels.update(1)
        print('\n')
        print('Design 1: Minimize Coolant Outlet Temperature')
        a = np.argmin(outlet_temp_arr_tot)
        print('Number of Channels: ' + str(num_channels_arr[a]))
        print('Coolant Outlet Temperature: ' + str(outlet_temp_arr_tot[a]))
        print('dP at min Outlet Temperature: ' + str(dP_arr_tot[a]))
        print(" ")
        print('Design 2: Minimize dP')
        a = np.argmin(dP_arr_tot)
        print('Number of Channels: ' + str(num_channels_arr[a]))
        print('Coolant Outlet Temperature: ' + str(outlet_temp_arr_tot[a]))
        print('dP at min Outlet Temperature: ' + str(dP_arr_tot[a]))

        plt.figure()
        plt.plot(num_channels_arr, dP_arr_tot)
        plt.xlabel('Number of Channels')
        plt.ylabel('dP')
        plt.figure()
        plt.xlabel('Number of Channels')
        plt.ylabel('Outlet Temp')
        plt.plot(num_channels_arr, outlet_temp_arr_tot)
        plt.show()
        return

    def simulate_regen(self):
        # Main loop

        if self.var_channel_w:
            self.var_width_optimization()
        else:
            self.fixed_width_optimization()

        return

    def single_design(self, num_channels, channel_w=.001):
        '''
        Calls either fixed_width_optimization() or var_width_optimization() for a single num channels and channel width (if fixed)
        Displays more specific plots of channel geometry, wall/coolant temps, and dP along nozzle contour
        Useful for debugging, analyzing trends, and pulling values for FEA/CFD, etc.
        '''
        if self.var_channel_w:
            self.num_channels = num_channels
            self.mDot_chan = self.engine.mDot_f / self.num_channels
            self.init_channel()

            x = 0
            while x < (self.numPTS - 1 - self.start_ind):
                i = self.numPTS - x - self.start_ind - 1  # Reverse order
                channel_w_high = ((self.engine.engineProps[
                                       i, 0] + self.wall_t) * np.pi * 2 - self.min_fin_w * self.num_channels) / self.num_channels
                channel_w_low = self.min_channel_w
                channel_w_range = np.arange(channel_w_low, channel_w_high, .0001)

                dP_arr = np.zeros(len(channel_w_range))
                dT_c_arr = np.zeros(len(channel_w_range))

                for k in np.arange(len(channel_w_range)):
                    (dP_arr[k], dT_c_arr[k]) = self.height_optimization(channel_w_range[k], i)
                a = np.min(dT_c_arr) * 1.02
                for l in np.arange(len(dT_c_arr)):
                    if dT_c_arr[l] > a:
                        dP_arr[l] = float('inf')
                a = np.argmin(dP_arr)

                self.height_optimization(channel_w_range[a], i, final=True)
                x += 1
            asdf = savgol_filter(self.channel_h_arr[:, 0], 31, 3)
            plt.figure()
            plt.plot(self.channel_h_arr)
            plt.plot(asdf)
            plt.show()
        else:

            self.num_channels = num_channels
            self.mDot_chan = self.engine.mDot_f / self.num_channels
            self.init_channel()
            x = 0
            while x < (self.numPTS - 1 - self.start_ind):
                i = self.numPTS - x - self.start_ind - 1  # Reverse order
                self.height_optimization(channel_w, i)
                x += 1


        plt.figure()
        plt.xlabel('Distance from Injector (m)')
        plt.ylabel('Dimension (m)')
        plt.title('Channel Geometry')
        plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.channel_h_arr)
        plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.channel_w_arr)
        plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.fin_w_arr)
        plt.legend(["Channel Height", "Channel Width", "Fin Width"])

        plt.figure()
        plt.title('Wall Temps')
        plt.xlabel('Distance from Injector (m)')
        plt.ylabel('Wall Temperature (K)')
        plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.wall_temps)
        plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.T_wc_arr)
        plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.coolant_temps)
        plt.legend(['Inner Wall Temp', 'Outer Wall Temp', 'Coolant Temp'])

        plt.figure()
        plt.title('Pressure Loss')
        plt.xlabel('Distance from Injector (m)')
        plt.ylabel('dP')
        plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.dP_arr)

        plt.show()

        return


'''
while np.abs(q1) > 1000:
    #This section solves for the required channel height to produce the required heat flux
    #Namely, solves q_out - q_in = 0 using the Newton-Raphson method
    channel_h_2 = channel_h_1 - .00001
    (h_c, A_effc, dP) = self.qout_calc(channel_h_2, channel_w_min, i)
    q2 = h_c * A_effc / length * (T_wc - self.T_cb) - q_in

    channel_h_1 = -q1 * (channel_h_1 - channel_h_2) / (q1 - q2) + channel_h_1
    if channel_h_1 <= 0:
        channel_h_1 = .0005
    (h_c, A_effc, dP) = self.qout_calc(channel_h_1, channel_w_min, i)
    q1 = h_c * A_effc / length * (T_wc - self.T_cb) - q_in
'''  # Newton-Raphson method (actually secant method)

