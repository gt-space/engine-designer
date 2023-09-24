import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from fuel_props import fuel_props
from bartz import bartz
from bartz import bartz_sigma
from scipy.optimize import fsolve
from scipy.signal import savgol_filter

from tqdm import tqdm
from structural import CopperWall
from contourDesigner.thermal_FEA import fea_solver

from contourDesigner.CEA_properties import get_props

from rocketcea.cea_obj import CEA_Obj

np.set_printoptions(threshold=sys.maxsize)  # Print full arrays (for debugging)

class RegenJacket:
    # Initialize the jacket object upon declaration
    def __init__(self, engine, wall_t=0.002, min_fin_w=0.001, min_channel_width=.001, var_channel_width=False,
                 T_ci=300,start_ind=0, min_channel_Re=4000, max_fin_AR = 5, max_channel_h = 0.01905, max_wall_temp = 1350):
        self.engine = engine  # Engine object to be jacketed

        self.wall_t = wall_t  # Inner wall thickness [m]
        self.wall_t_arr = np.ones([len(self.engine.engineProps[:, 0]) - 1, 1]) * self.wall_t
        self.min_fin_w = min_fin_w  # Minimum desired fin width [m]
        self.min_channel_w = min_channel_width

        self.var_channel_w = var_channel_width
        self.T_ci = T_ci  # Coolant inlet temperature [K]
        self.start_ind = start_ind  # Index where fuel enters (0 is exit plane)

        self.min_channel_Re = min_channel_Re
        self.max_fin_AR = max_fin_AR
        self.max_channel_h = max_channel_h

        self.numPTS = np.size(self.engine.engineProps, 0)  # Number of engine stations

        self.max_wall_temp_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1])

        self.cond_w = 375

        self.ax_pos = engine.engineProps[:-(self.start_ind + 1):, 1]

        self.max_wall_temp = max_wall_temp

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
        self.A_g_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1]) #Gas-side area array (length independent)
        self.A_effc_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1]) # Effective cooling area array (length independent)
        self.T_aw_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1]) # Adiabatic wall temp array
        self.length_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1]) # Length of station array

    def get_friction(self, Re, D_hyd):
        def colebrook(Re, D_hyd):

            k = 0.005 * (10 ** (-3))  # Surface roughness, guess based on our manufacuting capability, mm => m
            #Currently using Halbrook correlation to boost computational speed (generally accurate to ~3%)
            #To improve accuracy at significant computational expense, can solve numerically using bisection method
            '''
            fric_upper = 0.25
            fric_lower = 0.01
            fric_avg = 0.05
            fric_diff = -2*np.log10(k_OW/(3.71*D_hyd) + 2.51/(Re*np.sqrt(fric_avg))) - 1.0/np.sqrt(fric_avg)
            i = 0
            while np.abs(fric_diff)  > .00001 and i < 50:
                if fric_diff > 0:
                    fric_upper = fric_avg
                else:
                    fric_lower = fric_avg
                fric_avg = (fric_upper + fric_lower) / 2
                fric_diff = -2 * np.log10(k_OW / (3.71 * D_hyd) + 2.51 / (Re * np.sqrt(fric_avg))) - 1.0 / np.sqrt(fric_avg)
                i += 1


            fric = fric_avg
            if i == 50:
                fric = 1 / ((-1.8 * np.log10((k_OW / D_hyd / 3.7) ** 1.11 + 6.9 / Re)) ** 2)
            #func = lambda x: -2*np.log10(k_OW/(3.71*D_hyd) + 2.51/(Re*np.sqrt(x))) - 1.0/np.sqrt(x)
            #fric = fsolve(func, 1 / ((-1.8*np.log10((k_OW / D_hyd / 3.7)**1.11 + 6.9/Re))**2))
            '''


            fric = 1 / ((-1.8 * np.log10((k / D_hyd / 3.7) ** 1.11 + 6.9 / Re)) ** 2) #Initial guess: haaland equation

            cole_diff = -2*np.log10(k/D_hyd/3.7 + 2.51/Re/np.sqrt(fric)) - 1/np.sqrt(fric) #subtract LHS from RHS of colebrook
            while np.abs(cole_diff) > .00001:
                deriv = -2 /(np.log(10)*(k/D_hyd/3.7 + 2.51/Re/np.sqrt(fric))) * (2.51/Re*(-.5 * fric ** (-3/2))) + .5*fric**(-3/2)
                fric = fric - cole_diff/deriv #Newton Raphson Iteration
                cole_diff = -2 * np.log10(k / D_hyd / 3.7 + 2.51 / Re / np.sqrt(fric)) - 1 / np.sqrt(fric) #Recompute

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
        (rho_c, C_pc, cond_c, viscK_c) = fuel_props(self.T_cb)
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
        # viscK_c_w = fuel_props(T_wc)[3]  # Get viscosity at wall
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
        '''
        if i == 150 and channel_w > .001 and self.num_channels >= 60:
            ch_arr = np.arange(.0001, .01, .0001)
            eta_tot_arr = np.zeros(len(ch_arr))
            A_arr = np.zeros(len(ch_arr))
            for j in np.arange(len(ch_arr)):
                D_hyd = 2 * channel_w * ch_arr[j] / (channel_w + ch_arr[j])  # Channel hydraulic diameter
                A_cc = channel_w * ch_arr[j]  # Cross sectional area of channel (m^2)
                vel_c = self.mDot_chan / (A_cc * rho_c)  # Coolant velocity in channel (m/s)
                Re_c = (vel_c * D_hyd) / viscK_c  # Coolant Reynolds number
                Pr_c = viscK_c * rho_c * C_pc / cond_c  # Coolant Prandtl Number
                h_c = 0.023 * (Re_c ** 0.8) * (Pr_c ** 0.4) * (cond_c / D_hyd)  # Dittus Boelter reccomended for small dT Convection

                m = math.sqrt((2 * h_c) / (self.cond_w * fin_w))  # Fin efficiency length independent
                eta_fin = math.tanh(m * ch_arr[j]) / (m * ch_arr[j])  # Fin efficiency
                A_fin = 2 * ch_arr[j] * length  # Area of single fin in analysis segment
                A_finTot = self.num_channels * A_fin  # Total fin area
                A_wallTot = self.num_channels * channel_w * length  # Total wall contact area
                A_totc = A_finTot + A_wallTot  # Total coolant side heat transfer area
                eta_tot = (A_wallTot + A_finTot * eta_fin) / A_totc  # Overall convection efficiency

                A_effc = A_totc * eta_tot

                eta_tot_arr[j] = eta_tot
                A_arr[j] = A_effc
            plt.figure()

            plt.plot(ch_arr, eta_tot_arr)
            #plt.plot(ch_arr, A_arr)
            plt.xlabel('Channel Height, m')
            plt.ylabel('n_o')
            plt.title('Fin efficiency (.001m width, 60 channels)')
            plt.show()
            '''



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
        (rho_c, C_pc, cond_c, viscK_c) = fuel_props(self.T_cb)  # Get coolant fluid properties

        R = self.engine.engineProps[i, 0]  # Radius [m]
        length = self.engine.engineProps[i, 1] - self.engine.engineProps[i - 1, 1]  # Station Length [m]
        length = math.sqrt(length ** 2 + (R - self.engine.engineProps[i - 1, 0]) ** 2)  # Account for angle at station

        T_wg = self.max_wall_temp_arr[i - 1]  # DESIRED wall temp
        (h_g, q_in, T_aw) = bartz(self.engine, T_wg, i)  # Gas side heat flux parameters
        A_totg = 2 * math.pi * R * length  # total hot gas side wall area
        q_in = q_in * A_totg / length  # Heat flux to heat transfer (length independent bc some of the lengths are funky)
        T_wc = T_wg - (q_in / (self.cond_w / self.wall_t) / (A_totg / length))  # Fix later to add in variable wall cooling

        fin_w = ((2 * math.pi * (R + self.wall_t)) - (self.num_channels * channel_w)) / self.num_channels  # fin width


        channel_h_max_Re = 2 * self.mDot_chan / (self.min_channel_Re * rho_c * viscK_c) - channel_w  # Channel height upper bound, use as initial value for optimization
        channel_h_max_buckling = (2 * np.pi * R - self.num_channels * channel_w) / self.num_channels * self.max_fin_AR  # max fin height before buckling, as a multiple of fin width (FIX THIS LATER)


        channel_h_1 = min([channel_h_max_Re, channel_h_max_buckling, self.max_channel_h])  # Max channel height
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
            dP = float('inf')
            dT_c = float('inf')
            channel_h = channel_h_2
            q_out = q_in
            #return(float('inf'), float('inf'))

        elif q1 > 100: #q_out - q_in > 0 @ max channel height
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
            channel_h = channel_h_1
        else: #q_out = q_in somewhere within allowable channel range

            channel_h_avg = (channel_h_2 + channel_h_1) / 2
            (h_c, A_effc, dP) = self.qout_calc(channel_h_avg, channel_w, i)
            q_avg = h_c * A_effc / length * (T_wc - self.T_cb) - q_in
            while np.abs(q_avg) > 100:

                if q_avg > 0:
                    channel_h_2 = channel_h_avg
                else:
                    channel_h_1 = channel_h_avg
                channel_h_avg = (channel_h_1 + channel_h_2) / 2
                (h_c, A_effc, dP) = self.qout_calc(channel_h_avg, channel_w, i)
                q_avg = h_c * A_effc / length * (T_wc - self.T_cb) - q_in

            q1 = q_avg
            channel_h = channel_h_avg
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
            self.A_g_arr[i-1, 0] = 2 * np.pi * self.engine.engineProps[i-1, 0]
            self.A_effc_arr[i - 1, 0] = A_effc / length
            self.T_aw_arr[i - 1, 0] = T_aw
            self.length_arr[i - 1, 0] = length

        return (dP, dT_c)

    def fixed_width_optimization(self, num_channels_arr, channel_w_res):
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
        a = 2 * np.pi * np.min(self.engine.engineProps[:, 0]) / np.min(num_channels_arr)  # Filler variable for max channel width
        channel_w_range_max = np.arange(0, a + channel_w_res, channel_w_res)  # Array of overall possible channel widths to iterate through (e.g., for smallest num of channels)
        outlet_temp_3d = np.zeros([len(num_channels_arr), len(channel_w_range_max)]) #self.T_ci  # Coolant outlet temp for all possibilities #Edit 550 to make graph look pretty
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
                #while x < (self.numPTS - 1 - self.start_ind):
                for x in np.arange(self.numPTS - 1 - self.start_ind):
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
        outlet_temp_3d[outlet_temp_3d == 0] = np.min(outlet_temp_3d[outlet_temp_3d > 100]) - 10
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
        #plt.show()

        return (outlet_temp_3d, dP_tot_3d)

    def var_width_optimization(self, num_channels_arr):
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
        #plt.show()
        return

    def struct_optimzation(self, uts_FOS_min):
        x = 0
        T_cb = self.T_ci
        while x < (self.numPTS - 1 - self.start_ind):
            i = self.numPTS - x - self.start_ind - 1
            copper = CopperWall()
            T_wg_lower = 500
            (h_g, q_conv, T_aw) = bartz(self.engine, T_wg_lower, i)
            delT = q_conv / (self.cond_w / self.wall_t)
            uts_FOS_lower = copper.structural_analysis(delT, T_wg_lower)[2]

            T_wg_upper = 1250
            (h_g, q_conv, T_aw) = bartz(self.engine, T_wg_upper, i)
            delT = q_conv / (self.cond_w / self.wall_t)
            uts_FOS_upper = copper.structural_analysis(delT, T_wg_upper)[2]

            '''
            if i == self.engine.throatInd:
                temp_arr = np.arange(600,1300,10)
                uts_arr = np.zeros(len(temp_arr))
                yts_arr = np.zeros(len(temp_arr))
                sig_arr = np.zeros(len(temp_arr))
                for k_OW in np.arange(len(temp_arr)):
                    (h_g, q_conv, T_aw) = bartz(self.engine, temp_arr[k_OW], i)
                    delT = q_conv / (self.cond_w / self.wall_t)
                    (sig_equiv, yts_FOS, uts_FOS) = copper.structural_analysis(delT, temp_arr[k_OW])
                    uts_arr[k_OW] = uts_FOS
                    yts_arr[k_OW] = yts_FOS
                    sig_arr[k_OW] = sig_equiv
                plt.figure()
                plt.plot(temp_arr, uts_arr)
                plt.plot(temp_arr, yts_arr)
                plt.legend(["Ultimate", "Yield"])
                plt.xlabel("Wall Temperature (K)")
                plt.ylabel("Factor of Safety")
                plt.title("Wall Factor of Safety vs Temperature")


                plt.figure()
                plt.plot(temp_arr, sig_arr / 10**6)
                plt.plot(temp_arr, sig_arr * uts_arr / 10**6)
                plt.plot(temp_arr, sig_arr * yts_arr / 10**6)
                plt.legend(["Wall Stress", "Ultimate Stress", "Yield Stress"])
                plt.xlabel("Wall Temperature (K)")
                plt.ylabel("Stress (MPa)")

                plt.title("Wall Stresses vs Temperature")
                plt.show()
            '''
            '''T_wg = (T_wg_lower + T_wg_upper)/2'''
            T_wg = T_wg_lower + (uts_FOS_min - uts_FOS_lower)*(T_wg_upper - T_wg_lower)/(uts_FOS_upper - uts_FOS_lower)
            (h_g, q_conv, T_aw) = bartz(self.engine, T_wg, i)
            delT = q_conv / (self.cond_w / self.wall_t)
            uts_FOS = copper.structural_analysis(delT, T_wg)[2]

            while np.abs(uts_FOS - uts_FOS_min) > .001:
                if uts_FOS > uts_FOS_min:
                    T_wg_lower = T_wg
                    uts_FOS_lower = uts_FOS
                else:
                    T_wg_upper = T_wg
                    uts_FOS_upper = uts_FOS
                '''T_wg = (T_wg_lower + T_wg_upper) / 2'''
                T_wg = T_wg_lower + (uts_FOS_min - uts_FOS_lower)*(T_wg_upper - T_wg_lower)/(uts_FOS_upper - uts_FOS_lower)
                '''T_wg = (T_wg_lower*uts_FOS_upper - T_wg_upper*uts_FOS_lower)/(uts_FOS_upper - uts_FOS_lower)'''
                (h_g, q_conv, T_aw) = bartz(self.engine, T_wg, i)
                delT = q_conv / (self.cond_w / self.wall_t)
                (sig_equiv, yts_FOS, uts_FOS) = copper.structural_analysis(delT, T_wg)

            self.max_wall_temp_arr[i - 1] = np.min([T_wg, self.max_wall_temp])

            R = self.engine.engineProps[i, 0]
            length = self.engine.engineProps[i, 1] - self.engine.engineProps[i - 1, 1]  # Station Length [m]
            length = math.sqrt(
                length ** 2 + (R - self.engine.engineProps[i - 1, 0]) ** 2)  # Account for angle at station
            A_totg = 2 * np.pi * R * length
            q_conv = q_conv * A_totg
            (rho_c, C_pc, cond_c, viscK_c) = fuel_props(T_cb)
            T_cb += q_conv / self.engine.mDot_f / C_pc

            x += 1
        return T_cb

    def MR_optimization(self, outlet_temp, uts_FOS_min):
        print("Optimizing Mixture Ratio...")
        ispObj = CEA_Obj(oxName='LOX', fuelName='JetA', fac_CR=self.engine.con_rat)
        half = len(self.engine.engineContour[:, 0])
        half = int(half / 2)
        chBarrel = self.engine.engineContour[0:half, :]
        nozzleContour = self.engine.engineContour[half:self.engine.numPTS, :]
        MR_lower = 1
        MR_higher = 2.6
        MR = (MR_lower + MR_higher) / 2


        self.engine.engineProps = get_props(chBarrel, nozzleContour, self.engine.throatInd, ispObj, self.engine.P_inj_psi, MR, self.engine.A_t)
        T_cb = self.struct_optimzation(uts_FOS_min)
        i = 0
        while np.abs(T_cb - outlet_temp) > 0.1 and i < 100:
            if T_cb > outlet_temp:
                MR_higher = MR
            else:
                MR_lower = MR
            MR = (MR_lower + MR_higher) / 2
            self.engine.engineProps = get_props(chBarrel, nozzleContour, self.engine.throatInd, ispObj, self.engine.P_inj_psi, MR, self.engine.A_t)
            T_cb = self.struct_optimzation(uts_FOS_min)

            i += 1

        print("Required MR to close design: " + str(MR) + "\n")


        plt.figure()

        plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.max_wall_temp_arr)
        plt.xlabel("Distance from Injector (m)")
        plt.ylabel("Wall Temperature (K)")
        plt.title("Wall Temperatures to satisfy structural FOS")

        return MR

    def fea(self, i, res):
        #i = station
        corner_points = [(0,0),
                         (0, self.channel_h_arr[i,0] + self.wall_t),
                         (self.fin_w_arr[i,0]/2, self.channel_h_arr[i,0]+self.wall_t),
                         (self.fin_w_arr[i,0]/2, self.wall_t),
                         (self.fin_w_arr[i,0]/2+self.channel_w_arr[i,0]/2, self.wall_t),
                         (self.fin_w_arr[i,0]/2+self.channel_w_arr[i,0]/2, 0)]
        boundary_conds = np.array([
            [1, [corner_points[5], corner_points[0]], [self.h_g_arr[i], self.T_aw_arr[i]]],
            [1, [corner_points[2], corner_points[4]], [self.h_c_arr[i], self.coolant_temps[i]]],
        ], dtype=object)
        return fea_solver(corner_points, res, boundary_conds, self.cond_w)

    def simulate_regen(self, MR, num_channels_arr, channel_w_res):
        # Main loop
        yts_FOS_min = 0

        ispObj = CEA_Obj(oxName='LOX', fuelName='JetA', fac_CR=self.engine.con_rat)

        half = len(self.engine.engineContour[:, 0])
        half = int(half / 2)
        chBarrel = self.engine.engineContour[0:half, :]
        nozzleContour = self.engine.engineContour[half:self.engine.numPTS, :]
        self.engine.engineProps  = get_props(chBarrel, nozzleContour, self.engine.throatInd, ispObj, self.engine.P_inj_psi, MR, self.engine.A_t)



        if self.var_channel_w:
            self.var_width_optimization(num_channels_arr)
        else:
            self.fixed_width_optimization(num_channels_arr, channel_w_res)

        #self.film_cooling_2(100, .005)


        return

    def single_design(self, MR, num_channels, channel_w=.001, plot=False):
        '''
        Calls either fixed_width_optimization() or var_width_optimization() for a single num channels and channel width (if fixed)
        Displays more specific plots of channel geometry, wall/coolant temps, and dP along nozzle contour
        Useful for debugging, analyzing trends, and pulling values for FEA/CFD, etc.
        '''

        ispObj = CEA_Obj(oxName='LOX', fuelName='JetA', fac_CR=self.engine.con_rat)

        half = len(self.engine.engineContour[:, 0])
        half = int(half / 2)
        chBarrel = self.engine.engineContour[0:half, :]
        nozzleContour = self.engine.engineContour[half:self.engine.numPTS, :]
        self.engine.engineProps  = get_props(chBarrel, nozzleContour, self.engine.throatInd, ispObj, self.engine.P_inj_psi, MR, self.engine.A_t)

        if self.var_channel_w:
            self.num_channels = num_channels
            self.mDot_chan = self.engine.mDot_f / self.num_channels
            self.init_channel()

            x = 0
            while x < (self.numPTS - 1 - self.start_ind):
                print(x)
                i = self.numPTS - x - self.start_ind - 1  # Reverse order
                channel_w_high = ((self.engine.engineProps[
                                       i, 0] + self.wall_t) * np.pi * 2 - self.min_fin_w * self.num_channels) / self.num_channels
                channel_w_low = self.min_channel_w
                channel_w_range = np.arange(channel_w_low, channel_w_high, .0001)

                dP_arr = np.zeros(len(channel_w_range))
                dT_c_arr = np.zeros(len(channel_w_range))

                for k in np.arange(len(channel_w_range)):
                    (dP_arr[k], dT_c_arr[k]) = self.height_optimization(channel_w_range[k], i)
                '''
                a = np.min(dT_c_arr) * 1.05
                for l in np.arange(len(dT_c_arr)):
                    if dT_c_arr[l] > a:
                        dP_arr[l] = float('inf')
                '''
                a = np.argmin(dP_arr)

                self.height_optimization(channel_w_range[a], i, final=True)
                x += 1
            sf = savgol_filter(self.channel_h_arr[:, 0], 31, 3)
            plt.figure()
            plt.plot(self.channel_h_arr)
            plt.plot(sf)
            plt.xlabel("Station")
            plt.ylabel("Channel Height")
            plt.title("Channel Height for a variable-width channel")
            plt.legend(["Channel Height", "Curve-Smoothed"])
            #plt.show()

        else:

            self.num_channels = num_channels
            self.mDot_chan = self.engine.mDot_f / self.num_channels
            self.init_channel()
            x = 0
            while x < (self.numPTS - 1 - self.start_ind):
                i = self.numPTS - x - self.start_ind - 1  # Reverse order
                self.height_optimization(channel_w, i)
                x += 1

        copper = CopperWall()
        self.stress_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1])
        self.ult_FOS_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1])
        self.yield_FOS_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1])

        x = 0
        while x < (self.numPTS - 1 - self.start_ind):
            i = self.numPTS - x - self.start_ind - 1  # Reverse order
            T_wg = self.wall_temps[i-1]
            delT = T_wg - self.T_wc_arr[i-1]
            (stress, ult_FOS, yield_FOS) = copper.structural_analysis(delT, T_wg)
            self.stress_arr[i-1] = stress
            self.ult_FOS_arr[i-1] = ult_FOS
            self.yield_FOS_arr[i-1] = yield_FOS
            x += 1
        '''
        inds = [0, 100, self.engine.throatInd+100, 198]
        plt.figure()
        
        for i in range(len(inds)):
            area = (self.channel_h_arr[inds[i],0]+self.wall_t) * (self.channel_w_arr[inds[i],0]+self.fin_w_arr[inds[i],0])/2
            res = np.sqrt(area/10000)
            (mesh, mesh_points, temps) = self.fea(inds[i], res)

            x = np.array([mesh_points[j][0] for j in np.arange(len(mesh_points))])
            y = np.array([mesh_points[j][1] for j in np.arange(len(mesh_points))])

            temps = temps.astype('float')

            plt.subplot(2,2,i+1)
            plt.tricontourf(x, y, np.asarray(mesh.elements), temps, levels=100, cmap='jet')
            plt.colorbar()
            plt.axis('equal')
            title = ["Chamber Top", "Chamber Bottom", "Throat", "Nozzle Exit"]
            plt.title("Station %s: %s" % (inds[i], title[i]))
            plt.xlabel("Channel Width (m)")
            plt.ylabel("Channel Height (m)")
        '''
        '''
        plt.figure()
        plt.title('Wall Temps')
        plt.xlabel('Distance from Injector (m)')
        plt.ylabel('Wall Temperature (K)')
        plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.wall_temps)
        plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.T_wc_arr)
        plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.coolant_temps)
        plt.legend(['Inner Wall Temp', 'Outer Wall Temp', 'Coolant Temp'])
        '''


        if plot:
            plt.figure()
            plt.xlabel('Distance from Injector (m)')
            plt.ylabel('Dimension (m)')
            plt.title('Channel Geometry')
            plt.plot(self.ax_pos, self.channel_h_arr)
            plt.plot(self.ax_pos, self.channel_w_arr)
            plt.plot(self.ax_pos, self.fin_w_arr)
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
            plt.ylabel('dP (bar)')
            plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.dP_arr)

            plt.figure()
            plt.title('Factor of Safety')
            plt.xlabel('Distance from Injector (m)')
            plt.ylabel('FOS')

            plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.ult_FOS_arr)
            plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.yield_FOS_arr)
            plt.legend(["Ultimate", "Yield"])

            plt.figure()
            plt.title("Heat Flow into Wall")
            plt.xlabel("Distance from Injector (m)")
            plt.ylabel("Heat Flow (W)")
            plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.q_in_arr)
            plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.q_out_arr)

            plt.figure()
            plt.title("Heat Flux into Wall")
            plt.xlabel("Distance from Injector (m)")
            plt.ylabel("Heat Flow (W/m^2)")
            plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.q_in_arr/self.A_g_arr)
            plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.q_out_arr/self.A_effc_arr)
            #plt.show()

            plt.figure()
            plt.title("Heat Flux into Wall")
            plt.xlabel("Distance from Injector (m)")
            plt.ylabel("Heat Flow (W/m^2)")
            plt.plot(self.engine.engineProps[:-(self.start_ind + 1):, 1], self.h_g_arr)
            # plt.show()



        return

    def final_pass(self, MR, num_channels, c_h_arr, c_w_arr, wall_t):

        ispObj = CEA_Obj(oxName='LOX', fuelName='JetA', fac_CR=self.engine.con_rat)

        half = len(self.engine.engineContour[:, 0])
        half = int(half / 2)
        chBarrel = self.engine.engineContour[0:half, :]
        nozzleContour = self.engine.engineContour[half:self.engine.numPTS, :]
        self.engine.engineProps = get_props(chBarrel, nozzleContour, self.engine.throatInd, ispObj,
                                            self.engine.P_inj_psi, MR, self.engine.A_t)


        self.init_channel()
        self.num_channels = num_channels
        self.mDot_chan = self.engine.mDot_f / self.num_channels

        self.channel_h_arr = c_h_arr
        self.channel_w_arr = c_w_arr
        self.wall_t = wall_t

        copper = CopperWall()
        self.stress_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1])
        self.ult_FOS_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1])
        self.yield_FOS_arr = np.zeros([self.numPTS - 1 - self.start_ind, 1])

        x = 0

        while x < (self.numPTS - 1 - self.start_ind):
            i = self.numPTS - x - self.start_ind - 1

            (rho_c, C_pc, cond_c, viscK_c) = fuel_props(self.T_cb)  # Get coolant fluid properties

            R = self.engine.engineProps[i, 0]  # Radius [m]
            length = self.engine.engineProps[i, 1] - self.engine.engineProps[i - 1, 1]  # Station Length [m]
            length = math.sqrt(length ** 2 + (R - self.engine.engineProps[i - 1, 0]) ** 2)  # Account for angle at station

            A_totg = 2 * math.pi * R * length  # total hot gas side wall area
            fin_w = ((2 * math.pi * (R + self.wall_t)) - (self.num_channels * self.channel_w_arr[i-1])) / self.num_channels  # fin width


            (h_c, A_effc, dP) = self.qout_calc(self.channel_h_arr[i-1], self.channel_w_arr[i-1], i)

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

            dT_c = (q_in + q_out) * length / 2 / (C_pc * self.engine.mDot_f)  # Temperature change in coolant


            self.T_cb += dT_c
            self.P_co -= dP

            self.fin_w_arr[i - 1, 0] = fin_w
            self.wall_temps[i - 1, 0] = T_wg
            self.coolant_temps[i - 1, 0] = self.T_cb
            self.pressures[i - 1, 0] = self.P_co
            self.h_g_arr[i - 1, 0] = h_g
            self.h_c_arr[i - 1, 0] = h_c
            self.T_wc_arr[i - 1, 0] = T_wc
            self.dP_arr[i - 1, 0] = dP
            self.Re_c_arr[i - 1, 0] = self.mDot_chan / (rho_c * viscK_c) * 2 / (
                    self.channel_w_arr[i-1] + self.channel_h_arr[i-1])  # Double check that this is right
            self.q_in_arr[i - 1, 0] = q_in
            self.q_out_arr[i - 1, 0] = q_out
            self.A_effc_arr[i - 1, 0] = A_effc / length
            self.T_aw_arr[i - 1, 0] = T_aw
            self.length_arr[i - 1, 0] = length


            delT = T_wg - self.T_wc_arr[i-1]
            (stress, yield_FOS, ult_FOS) = copper.structural_analysis(delT, T_wg)

            self.stress_arr[i-1] = stress
            self.ult_FOS_arr[i-1] = ult_FOS
            self.yield_FOS_arr[i-1] = yield_FOS

            x += 1

        return

'''
    def struct_maximization(self):
        phi = (1+np.sqrt(5))/2
        x = 0
        T_cb = self.T_ci
        while x < (self.numPTS - 1 - self.start_ind):
            i = self.numPTS - x - self.start_ind - 1
            copper = CopperWall()
            T_wg_lower = 500
            (h_g, q_conv, T_aw) = bartz(self.engine, T_wg_lower, i)
            delT = q_conv / (self.cond_w / self.wall_t)
            uts_FOS_lower = copper.structural_analysis(delT, T_wg_lower)[2]

            T_wg_upper = 1500
            (h_g, q_conv, T_aw) = bartz(self.engine, T_wg_upper, i)
            delT = q_conv / (self.cond_w / self.wall_t)
            uts_FOS_upper = copper.structural_analysis(delT, T_wg_upper)[2]

            T_wg_midlow = T_wg_lower + (T_wg_upper-T_wg_lower)*phi/(2*phi+1)
            (h_g, q_conv, T_aw) = bartz(self.engine, T_wg_midlow, i)
            delT = q_conv / (self.cond_w / self.wall_t)
            uts_FOS_midlow = copper.structural_analysis(delT, T_wg_midlow)[2]

            T_wg_midup = T_wg_lower + (T_wg_upper-T_wg_lower)*(phi+1)/(2*phi+1)
            (h_g, q_conv, T_aw) = bartz(self.engine, T_wg_midup, i)
            delT = q_conv / (self.cond_w / self.wall_t)
            uts_FOS_midup = copper.structural_analysis(delT, T_wg_midup)[2]


            if uts_FOS_upper > uts_FOS_midup:
                T_wg = T_wg_upper
            elif uts_FOS_lower > uts_FOS_midlow:
                T_wg = T_wg_lower
            else:

                while np.abs(uts_FOS_midup - uts_FOS_midlow) > 0.000001:
                    if uts_FOS_midup > uts_FOS_midlow:
                        T_wg_lower = T_wg_midlow; uts_FOS_lower = uts_FOS_midlow
                        T_wg_midlow = T_wg_midup; uts_FOW_midlow = uts_FOS_midup

                        T_wg_midup = T_wg_lower + (T_wg_upper-T_wg_lower)*(phi+1)/(2*phi+1)
                        (h_g, q_conv, T_aw) = bartz(self.engine, T_wg_midup, i)
                        delT = q_conv / (self.cond_w / self.wall_t)
                        uts_FOS_midup = copper.structural_analysis(delT, T_wg_midup)[2]
                    elif uts_FOS_midup < uts_FOS_midlow:
                        T_wg_upper = T_wg_midup; uts_FOS_upper = uts_FOS_midup
                        T_wg_midup = T_wg_midlow; uts_FOS_midup = uts_FOS_midlow

                        T_wg_midlow = T_wg_lower + (T_wg_upper-T_wg_lower)*phi/(2*phi+1)
                        (h_g, q_conv, T_aw) = bartz(self.engine, T_wg_midlow, i)
                        delT = q_conv / (self.cond_w / self.wall_t)
                        uts_FOS_midup = copper.structural_analysis(delT, T_wg_midlow)[2]
                T_wg = (T_wg_midlow+T_wg_midup)/2

            self.max_wall_temp_arr[i - 1] = T_wg

            R = self.engine.engineProps[i, 0]
            length = self.engine.engineProps[i, 1] - self.engine.engineProps[i - 1, 1]  # Station Length [m]
            length = math.sqrt(
                length ** 2 + (R - self.engine.engineProps[i - 1, 0]) ** 2)  # Account for angle at station
            A_totg = 2 * np.pi * R * length
            q_conv = q_conv * A_totg
            (rho_c, C_pc, cond_c, viscK_c) = fuel_props(T_cb)
            T_cb += q_conv / self.engine.mDot_f / C_pc

            x += 1
        return
'''