import numpy as np
from matplotlib import pyplot as plt

from bartz import bartz

class HeatSink:
    def __init__(self, engine, run_time, T_wall_i, split):
        self.rad_c = 2.2 * .0254
        self.D_c = engine.engineProps[0, 0] * 2
        self.D_t = np.min(engine.engineProps[:, 0]) * 2
        self.P_0 = engine.engineProps[0, 8] * 100000
        self.mdot = engine.mDot_tot
        self.c_star = engine.C_star
        self.run_time = run_time
        self.T_wall_i = T_wall_i #K, initial temp
        self.split = split
        self.engine = engine

    class Steel:
        def __init__(self):
            self.k = 42.7 # W/m*K
            self.rho = 7850 # kg/m3
            self.Cp = 477 # J/kg*K
            #self.ind = np.arange(split)

    class Graphite:
        def __init__(self):
            self.k = 130  # W/m*K
            self.rho = 1760  # kg/m3
            self.Cp = 700 # J/kg*K
            #self.ind = np.arange(split, len(engine.engineProps[:, 0]))
            #siphen:::

            self.k = .225  # W/m*K
            self.rho = 1716  # kg/m3
            self.Cp = 1000  # J/kg*K
            self.CTE = 1.368e-5  # m/m/C = m/m/K
            self.v = .3575  # Google
            self.E = 22e9  # Sutton
            self.UTS = 6.2e7  # Pa Matweb

            '''class SilicaPhenolic:
                def __init__(self):
                    self.k_OW = .225  # W/m*K
                    self.rho_OW = 1716  # kg/m3
                    self.Cp_OW = 1000  # J/kg*K
                    self.CTE_OW = 1.368e-5  # m/m/C = m/m/K
                    self.v = .3575  # Google
                    self.E = 22e9  # Sutton
                    self.UTS_OW = 6.2e7  # Pa Matweb'''



    class Copper:
        def __init__(self):
            self.k = 350 #W/mK
            self.rho = 8000
            self.Cp = 385
            #self.ind = np.arange(split)

    def steel_bounds(self, T_arr, temp_lim):
        """Determine points where converging section steel must end and diverging section steel can start"""
        i = 0
        while T_arr[i] <= temp_lim: i += 1
        self.conv_end = i

        j = len(self.engine.engineProps[:,0]) - 1
        while T_arr[j] <= temp_lim: j += -1
        self.div_start = j
        return i, j



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

        zeta_arr = zeta_solver(15, Bi)
        Cn_arr = 4 * np.sin(zeta_arr) / (2*zeta_arr + np.sin(2*zeta_arr))

        theta_0 = np.zeros(len(Fo_arr)) #Theta = T(t) - T_wg / T_wall_i - T_wg, T(t) is OUTER wall temp
        theta_end = np.zeros(len(Fo_arr)) #Theta = T(t) - T_wg / T_wall_i - T_wg, T(t) is INNER wall temp
        for i in np.arange(len(Fo_arr)):
            for j in np.arange(len(zeta_arr)):
                theta_0[i] += Cn_arr[j] * np.exp(-1 * zeta_arr[j] **2 * Fo_arr[i])
                theta_end[i] += Cn_arr[j] * np.exp(-1 * zeta_arr[j] **2 * Fo_arr[i]) * np.cos(zeta_arr[j] * x_star)
        T_wall_inner = theta_end * (T_wall - T_ad) + T_ad #Inner wall temp - this is what matters

        return T_wall_inner

    def chamber_analysis(self, run_time, iter = 5, plot2d = False, plot3d = False):
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
                plt.plot(self.engine.engineContour[:,1], T_arr_3d[-1, :])
                plt.xlabel('Axial Station')
                plt.ylabel('Wall Temp (K)')
                plt.title('Wall Temp vs. Axial Station')

                plt.figure()
                plt.plot(self.engine.engineContour[:,1], self.engine.engineProps[:, 9])
                #plt.show()
            print("Chamber Wall Temp, " + str(run_time) + " seconds: ", T_arr_3d[-1, 0])
        #Steel section:

        steel = self.Steel()
        gr = self.Graphite()
        time_arr = np.arange(.001, run_time+.01, 0.01)
        T_arr_3d = np.empty([len(time_arr), len(self.engine.engineProps[:, 0])], float)
        h_g = np.zeros(len(self.engine.engineProps[:, 0]))
        for i in np.arange(len(self.engine.engineProps[:, 0])):
            if i < self.split:
                mat = steel
                thickness = 1.5 * .0254
            else:
                mat = gr
                thickness = (0.4 * .0254) + self.engine.engineProps[self.split, 0] - self.engine.engineProps[i, 0]
            h_g[i] = bartz(self.engine, self.T_wall_i, i)[0]

            Fo_arr = mat.k / mat.rho / mat.Cp * time_arr / thickness ** 2
            for z in np.arange(iter):
                Bi = h_g[i] * thickness / mat.k
                T_arr = self.transient_wall_temp(Bi, Fo_arr, self.T_wall_i, self.engine.engineProps[i, 9])
                h_g[i] = bartz(self.engine, np.mean(T_arr), i)[0]
            T_arr_3d[:, i] = T_arr


        if plot2d or plot3d:
            plot_temp()
        return T_arr_3d

    def wall_temp_gradient(self, station, thickness, mat):
        x_star_arr = np.arange(0, 1.01, .01)
        x_star_arr = np.flip(x_star_arr)
        time_arr = np.arange(.01, self.run_time + .01, 0.01)
        h_g = bartz(self.engine, self.T_wall_i, station)[0]
        Fo_arr = mat.k_OW / mat.rho_OW / mat.Cp_OW * time_arr / thickness ** 2
        T_arr = np.zeros(len(x_star_arr))
        for i in np.arange(len(x_star_arr)):
            if i == 0:
                for z in np.arange(5):
                    Bi = h_g * thickness / mat.k_OW
                    T = self.transient_wall_temp(Bi, Fo_arr, self.T_wall_i, self.engine.engineProps[station, 9], x_star_arr[i])
                    h_g = bartz(self.engine, np.mean(T) , station)[0]
            else:
                Bi = h_g * thickness / mat.k_OW
                T = self.transient_wall_temp(Bi, Fo_arr, self.T_wall_i, self.engine.engineProps[station, 9], x_star_arr[i])
            T_arr[i] = T[-1]
        T_arr = np.flip(T_arr)

        plt.figure()
        plt.plot(x_star_arr * thickness / .0254, T_arr)
        plt.xlabel('Thickness (in)')
        plt.ylabel('Wall Temp (K)')
        plt.title('Temperature across Wall Thickness, Chamber')
        #plt.show()





'''
rad_arr = np.linspace(hs.rad[station], hs.rad[station]+thickness, len(T_arr))
u = 12.2*10**(-6)
rad_arr_new = rad_arr * (T_arr - 300) * u + rad_arr

plt.figure()
plt.plot(rad_arr)
plt.plot(rad_arr_new)
plt.title('Radius Correction for Thermal Expansion')
plt.ylabel('Radial Distance from Centerline (in)')
#plt.show()
'''
'''
wall_temp_gradient(heatsink, 100, 1 * .0254, heatsink.Steel())


def conv_coeff_vs_wall_temp(hs, station):
    T_wall_i = 300 #K
    T_wall_max = 1500 #K
    T_arr = np.arange(T_wall_i, T_wall_max, 5)
    conv_coeff_arr = np.empty(len(T_arr), float)
    for i in np.arange(len(conv_coeff_arr)):
        conv_coeff_arr[i] = bartz(engine, T_arr[i], station)[0]

    plt.figure()
    plt.plot(T_arr, conv_coeff_arr)
    plt.xlabel('Wall Temp (K)')
    plt.ylabel('Conv. Coeff (W/m2K)')
    plt.title(' Convection Coeff. vs Wall Temp')

    #plt.show()

conv_coeff_vs_wall_temp(heatsink, 1)
'''
plt.show()


