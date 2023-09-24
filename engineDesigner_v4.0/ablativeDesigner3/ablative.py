import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve
from design import Engine
from bartz import bartz
import icecream as ic
# This script designs the self engine. It takes in an engine class, plus a run time, and outputs an engine
# contour, liner thickness, and overwrap thickness

class Ablative:
    def __init__(self, engine, run_time, T_wall_i, split):
        self.engine = engine

        # Chamber Properties
        self.D_c = engine.engineProps[0, 0] * 2
        self.P_0 = engine.engineProps[0, 8] * 100000  # Pa
        self.T_0 = engine.engineProps[0, 9]  # K
        self.mdot = engine.mDot_tot
        self.c_star = engine.C_star
        self.D_t = np.min(engine.engineProps[:, 0]) * 2
        self.rad_c = 2.2 * .0254

        # Exit Properties
        self.r_e = engine.engineProps[-1, 0]
        self.MW_e = engine.engineProps[-1, 13]
        self.gamma_e = engine.engineProps[-1, 15]

        # Ablation Properties
        self.run_time = run_time

        # Thermal Transient Property data
        self.T_wall_i = T_wall_i #K, initial temp
        self.split = split


    class Steel:
        def __init__(self):
            self.k = 42.7 # W/m*K
            self.rho = 7850 # kg/m3
            self.Cp = 477 # J/kg*K
            #self.ind = np.arange(split)
    class SilicaPhenolic:
        def __init__(self):
            self.k = .225 # W/m*K
            self.rho = 1716 # kg/m3
            self.Cp = 1000 # J/kg*K
            self.CTE = 1.368e-5 # m/m/C = m/m/K
            self.v = .3575 # Google
            self.E = 22e9 # Sutton
            self.UTS = 6.2e7 #Pa Matweb

    class SilicaEpoxy:  # matweb
        def __init__(self):
            self.k = .225  # W/m*K #number from literature
            self.rho = 2300  #
            self.Cp = 1000  #
            self.CTE = 17e-6  #
            self.v = .23  #
            self.E = 36e9  #
            self.UTS = 15.2e6  #
            #self.ind = np.arange(split)
    class CarbonEpoxy:
        def __init__(self):
            self.k = 10 # conservative thing https://www.researchgate.net/profile/Ronald-Joven/publication/288102626_Thermal_properties_of_carbon_fiber-epoxy_composites_with_different_fabric_weaves/links/56be27e408aee5caccf2f5d3/Thermal-properties-of-carbon-fiber-epoxy-composites-with-different-fabric-weaves.pdf
            self.rho = 1410  # kg/m3
            self.Cp = 1200 # J/kg*K
            self.v_OW = .4  # using epoxy value
            self.E_OW = 7e10  # Pa Matweb
            self.CTE = 2.13e-5 #m/m/C = m/m/K
            self.v = .4
            self.E = 7e10
            self.UTS = 918e6 #Pa Matweb

            #self.ind = np.arange(split)
    class Graphite:
        def __init__(self):
            self.k = 110  # W/m*K
            self.rho = 1760  # kg/m3
            self.Cp = 700 # J/kg*K
            self.UTS = 4200*6894.76 # Pa really i s just flexural strength in tension
            #self.ind = np.arange(split, len(engine.engineProps[:, 0]))

    class Copper:
        def __init__(self):
            self.k = 350 #W/mK
            self.rho = 8000
            self.Cp = 385
            #self.ind = np.arange(split)

    def closest(self,lst, K):
        lst = np.asarray(lst)
        idx = (np.abs(lst - K)).argmin()
        return idx

    def AonAStar_value(self, M, gamma):  # value of AonAstar for a given Mach and gamma
        return (1 / M) * ((1 + ((gamma - 1) / 2) * (M ** 2)) / ((gamma + 1) / 2)) ** (
                    (gamma + 1) / (2 * (gamma - 1)))

    def AonAStar_mach(self, Arat, gamma):  # mach number for a given area ratio
        def A_over_A_star(M):
            return (1 / M) * ((1 + ((gamma - 1) / 2) * (M ** 2)) / ((gamma + 1) / 2)) ** (
                        (gamma + 1) / (2 * (gamma - 1))) - Arat

        subsonic_M = fsolve(A_over_A_star, 0.001)
        supersonic_M = fsolve(A_over_A_star, 20)
        # print('sub: ', subsonic_M, '\n super: ', supersonic_M)
        return subsonic_M, supersonic_M


    def compflowtool(self,gamma, mach):

        p_0_over_p = (1 + ((gamma - 1) / 2) * mach ** 2) ** (gamma / (gamma - 1))
        T_0_over_T = (1 + ((gamma - 1) / 2) * mach ** 2)
        rho_0_over_rho = (1 + ((gamma - 1) / 2) * mach ** 2) ** (1 / (gamma - 1))
        A_over_Astar = (1 / mach) * ((1 + ((gamma - 1) / 2) * mach ** 2) / ((gamma + 1) / 2)) ** (
                    (gamma + 1) / (2 * (gamma - 1)))
        MFP = (np.sqrt(gamma) * mach) / (1 + ((gamma - 1) / 2) * mach ** 2) ** ((gamma + 1) / (2 * (gamma - 1)))

        return p_0_over_p, T_0_over_T, rho_0_over_rho, A_over_Astar, MFP


    def steel_bounds(self, T_arr, temp_lim):
        """Determine points where converging section steel must end and diverging section steel can start"""
        i = 0
        while T_arr[i] <= temp_lim: i += 1
        self.conv_end = i

        j = len(self.engine.engineProps[:,0]) - 1
        while T_arr[j] <= temp_lim: j += -1
        self.div_start = j
        return i, j


# ablative wall thickness sizing
    def size_ablator_thickness(self, burn_time,T_w, T_amb=272,wall_idx=[],t_A_design=[],liner=[],OW=[],I=[]):
        # Size ablative thickness if not already specified:
        if t_A_design == []:
            # First order empirical correlation from NASA SP8124:
            t_A_NASA =0.04 *burn_time**(1/2) * .0254
            # Using Sutton "Typical values" for erosion rate
            t_A_Sutton = .015*burn_time*.0254
            t_A = (t_A_NASA+t_A_Sutton)/2
        else:
            t_A = t_A_design

        # Tamb is # 30F since colder is more conservative
        if wall_idx == []:
            wall_idx = self.engine.throatInd


        if OW ==[]:
            OW= self.CarbonEpoxy()
        if I ==[]:
            I= self.SilicaEpoxy()
        if liner==[]:
            liner= self.SilicaPhenolic()

        #Material properties:
        # Ablator:
        k_A = liner.k  # W/mK

        # Structural overwrap (OW):
        k_OW = OW.k  # W/mK Matweb

        (h_g, q_conv, T_aw) = bartz(self.engine, T_w, wall_idx)
        # Difference between T_w and TadW?

        T_I = T_w - q_conv * t_A / k_A
        h_air = 100  # higher is more conservative

        T_OW = q_conv / h_air + T_amb
        #T_OW = 273.15
        t_OW = k_OW * (T_I - T_OW) / q_conv

        # Wall Temperature Plots

        return t_A,t_OW,q_conv,h_g,T_aw,T_I,T_OW


    def estim_liner_mass(self,engine,t_w,t_I,t_OW,mat_w,mat_I,mat_OW):# t_w, t_I, t_OW are in METERS!
        x = engine.engineProps[:, 1]
        y = engine.engineProps[:, 0]


        def compute_hollow_cyl_vol(x, y, R, r): # compute
            # Calculate the distance between each cylindrical cell
            cell_dist = np.sqrt((x[1:] - x[:-1]) ** 2 + (y[1:] - y[:-1]) ** 2)

            # Calculate the volume of each cylindrical element
            shell_volumes = np.pi * ((R[:-1]) ** 2 - (r[:-1]) ** 2) * cell_dist

            # Compute the total layer volume by summing up the volumes of all shells
            layer_volume = np.sum(shell_volumes)
            return layer_volume

        # Calculate the outer and inner radii of the cylinder
        # For the LINER: (W)
        R = y + t_w
        r = y
        vol_w = compute_hollow_cyl_vol(x, y, R, r)
        # for the INSULATION (I)
        R = y + t_w + t_I
        r = y + t_w
        vol_I = compute_hollow_cyl_vol(x, y, R, r)
        # for the OVERWRAP (OW)
        R = y + t_w + t_I + t_OW
        r = y + t_w + t_I
        vol_OW = compute_hollow_cyl_vol(x, y, R, r)
        mass = vol_w * mat_w.rho_OW + vol_I * mat_I.rho_OW + vol_OW * mat_OW.rho_OW
        mass_lbm = 2.20462*mass
        return mass_lbm


    def wall_temp_plot(self,T_w,mat,t,T_amb=272,numpts=200,wall_idx=[]):

        if wall_idx == []:
            wall_idx = self.engine.throatInd
        (h_g, q_conv, T_aw) = bartz(self.engine, T_w, wall_idx)


        k = mat.k_OW


    def chamber_stress_analysis(self,material,t,DT):
        p_c = self.engine.P_inj*100000
        r_c = self.engine.engineProps[0, 0]
        CTE = material.CTE_OW
        E = material.E
        v = material.v
        UTS=material.UTS_OW

        stress_Hoop = p_c*r_c/t #Pa
        stress_Thermal =2*CTE*E*DT/(1-v)
        stress = stress_Hoop+stress_Thermal
        FOS = UTS/stress
        return stress,FOS
    def throat_insert_size(self,material,FOS):
        #choosing chamber radius and pressure at end of converging section
        p_max =self.P_0 #Pa
        r_max =self.engine.engineContour[0,0] #m
        t_th= p_max*r_max/(material.UTS_OW / FOS) #meters
        return t_th # in METERS


    def transient_wall_temp(self, Bi, Fo_arr, T_wall, T_ad, x_star = 1):
        """Solves for the outer and inner wall temperatures for a given run time and material properties"""
        def zeta_solver(num, Bi):
            """Solves for the first n roots of the equation zeta * tan(zeta) = Bi"""
            #Midpoint method for solving for zeta
            zeta_arr = np.zeros(num)
            for n in np.arange(num):
                high = np.pi * n + np.pi/2
                low = np.pi * n
                avg = (high + low) /2
                err = avg * np.tan(avg) - Bi
                while np.abs(err) >= .0001:
                    if err > 0:
                        high = avg
                    else:
                        low = avg
                    avg = (high + low) / 2
                    err = avg * np.tan(avg) - Bi
                zeta_arr[n] = np.abs(avg)
            return zeta_arr

        zeta_arr = zeta_solver(7, Bi)
        Cn_arr = 4 * np.sin(zeta_arr) / (2*zeta_arr + np.sin(2*zeta_arr))

        theta_0 = np.zeros(len(Fo_arr)) #Theta = T(t) - T_wg / T_wall_i - T_wg, T(t) is OUTER wall temp
        theta_end = np.zeros(len(Fo_arr)) #Theta = T(t) - T_wg / T_wall_i - T_wg, T(t) is INNER wall temp
        for i in np.arange(len(Fo_arr)):
            for j in np.arange(len(zeta_arr)):
                theta_0[i] += Cn_arr[j] * np.exp(-1 * zeta_arr[j] **2 * Fo_arr[i])
                theta_end[i] += Cn_arr[j] * np.exp(-1 * zeta_arr[j] **2 * Fo_arr[i]) * np.cos(zeta_arr[j] * x_star)
        T_wall_inner = theta_end * (T_wall - T_ad) + T_ad #Inner wall temp - this is what matters

        for i in np.arange(len(Fo_arr)):
            for j in np.arange(len(zeta_arr)):
                theta_0[i] += Cn_arr[j] * np.exp(-1 * zeta_arr[j] **2 * Fo_arr[i])
                theta_0[i] += Cn_arr[j] * np.exp(-1 * zeta_arr[j] **2 * Fo_arr[i]) * np.cos(zeta_arr[j] * x_star)
        T_wall_outer = theta_0 * (T_wall - T_ad) + T_ad #Inner wall temp - this is what matters


        return T_wall_inner,T_wall_outer

    def chamber_analysis(self, run_time,liner_thickness, iter = 5, plot2d = False, plot3d = False,T='inner'):
        if T == 'outer':  # display temp gradient for outer wall temp
            T = 1
        else:
            T = 0

        def plot_temp():
            if plot3d:
                X, Y = np.meshgrid(np.arange(len(self.engine.engineProps[:, 0])), time_arr)
                # Plot data
                plt.figure()
                ax_ch = plt.axes(projection='3d')
                ax_ch.plot_surface(X, Y, T_arr_3d, rstride=1, cstride=1,
                               cmap='viridis', edgecolor='none')
                ax_ch.set_xlabel('Axial Station')
                ax_ch.set_ylabel('Time (sec)')
                ax_ch.set_zlabel('Wall Temp (K)')
                plt.title('Wall Temperature over time')
                #plt.show()
            if plot2d:
                plt.figure()
                plt.plot(self.engine.engineContour[:,1]/.0254, T_arr_3d[-1, :])
                plt.xlabel('Axial Station dist (in)')
                plt.ylabel('Wall Temp (K)')
                plt.title('Wall Temp vs. Axial Station')

                #plt.figure()
                #plt.plot(self.engine.engineContour[:,1]/.0254, self.engine.engineProps[:, 9])
                #plt.xlabel('Axial Station dist (in)')
                #plt.ylabel('Pressure (')
                #plt.title('Wall Temp vs. Axial Station')
                #plt.show()
            print("Chamber Wall Temp, " + str(run_time) + " seconds: ", T_arr_3d[-1, 0])
        #CF OW section:
        OW = self.CarbonEpoxy()
        liner = self.SilicaPhenolic()
        time_arr = np.arange(.001, run_time+.01, 0.1)
        T_arr_3d = np.empty([len(time_arr), len(self.engine.engineProps[:, 0])], float)
        h_g = np.zeros(len(self.engine.engineProps[:, 0]))
        for i in np.arange(len(self.engine.engineProps[:, 0])):
            if i < self.split:
                mat = OW
                thickness = 1.5 * .0254
            else:
                mat = liner
                thickness = (liner_thickness) + self.engine.engineProps[self.split, 0] - self.engine.engineProps[i, 0]
            h_g[i] = bartz(self.engine, self.T_wall_i, i)[0]

            Fo_arr = mat.k / mat.rho / mat.Cp * time_arr / thickness ** 2
            for z in np.arange(iter):
                Bi = h_g[i] * thickness / mat.k
                T_arr = self.transient_wall_temp(Bi, Fo_arr, self.T_wall_i, self.engine.engineProps[i, 9])[T]
                h_g[i] = bartz(self.engine, np.mean(T_arr), i)[0]
            T_arr_3d[:, i] = T_arr


        if plot2d or plot3d:
            plot_temp()
        return T_arr_3d

    def wall_temp_gradient(self, station, thickness, mat,T='inner'):
        if T=='outer': # display temp gradient for outer wall temp
            T=1
        else:
            T=0
        x_star_arr = np.arange(0, 1.01, .01)
        x_star_arr = np.flip(x_star_arr)
        time_arr = np.arange(.01, self.run_time + .01, 0.01)
        h_g = bartz(self.engine, self.T_wall_i, station)[0]
        Fo_arr = mat.k_OW / mat.rho_OW / mat.Cp_OW * time_arr / thickness ** 2
        T_arr = np.zeros(len(x_star_arr))
        for i in np.arange(len(x_star_arr)):
            if i == 0:
                for z in np.arange(3):
                    Bi = h_g * thickness / mat.k_OW
                    T = self.transient_wall_temp(Bi, Fo_arr, self.T_wall_i, self.engine.engineProps[station, 9], x_star_arr[i])[T]
                    h_g = bartz(self.engine, np.mean(T) , station)[0]
            else:
                Bi = h_g * thickness / mat.k_OW
                T = self.transient_wall_temp(Bi, Fo_arr, self.T_wall_i, self.engine.engineProps[station, 9], x_star_arr[i])[T]
            T_arr[i] = T[-1]
        T_arr = np.flip(T_arr)

        plt.figure()
        plt.plot(x_star_arr * thickness / .0254, T_arr)
        plt.xlabel('Thickness (in)')
        plt.ylabel('Wall Temp (K)')
        plt.title('Temperature across Wall Thickness, Chamber')
        plt.show()


plt.show()


