import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve
from design import Engine
from bartz import bartz
import icecream as ic

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

    class SilicaEpoxy: #matweb
        def __init__(self):
            self.k = .225 # W/m*K #number from literature
            self.rho = 2300  #
            self.Cp = 1000  #
            self.CTE = 17e-6  #
            self.v = .23  #
            self.E = 36e9 #
            self.UTS = 15.2e6  #

            #self.ind = np.arange(split)
    class CarbonEpoxy:
        def __init__(self):
            self.k = .8 # W/m*K #using middle range (around 150C) varies between 800 and 1400
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

    def thrust_curve(self, throat_insert_ablation_rate,chamber_ablation_rate, burnTime,debugMach,showplot=True):

        burnTimes = np.arange(0,burnTime,.01)
        Thrust = np.zeros(burnTimes.size)
        Thrust_mdot_cnst = np.zeros(burnTimes.size)
        M_e= np.zeros(burnTimes.size)
        for i in range(len(burnTimes)):
            (Thrust_mdot_cnst[i],Thrust[i],M_e[i]) = self.thrust_change_throat(throat_insert_ablation_rate,chamber_ablation_rate, burnTimes[i])
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
            plt.plot(burnTimes, Thrust, 'b*',label='Thrust with changing Mdot')
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
        r_t = self.D_t/2
        eroded_r_t = r_t - Abl_rate_thr*time
        eroded_A_t = np.pi*eroded_r_t**2
        eroded_epsilon = A_e/eroded_A_t

        # Eroded engine properties

        M_e = self.AonAStar_mach(eroded_epsilon, self.gamma_e)[1]

        p_0_over_p, T_0_over_T, rho_0_over_rho, A_over_Astar, MFP = self.compflowtool(self.gamma_e, M_e)
        T_e = T_0_over_T ** -1 * self.T_0
        R_e = 8314.5/self.MW_e
        u_e = M_e * np.sqrt(self.gamma_e * R_e * T_e)

        # if mass flow rate doesnt change:
        Thrust_mdot_cnst = self.mdot*u_e

        # if mass flow rate DOES change as a result of chamber pressure dropping:


        gamma_c = self.engine.engineProps[0,15]
        gamma = (gamma_c+self.gamma_e)/2
        MW_c = self.engine.engineProps[0,13]
        MW = (MW_c+self.MW_e)/2
        R = 8314.5/(MW)
        p_c_new = self.p_c_drop(chamber_ablation_rate,time,showplot=False)
        A = p_c_new/np.sqrt(R*self.T_0) *eroded_A_t
        B = np.sqrt(gamma*(2/(gamma+1))**((gamma+1)/(gamma-1)))
        newmdot = A * B
        Thrust = newmdot * u_e
        # This assumes the following:
        # Chamber temperature does not change due to ablation
        # Molecular Weight (MW) and specific heat ratio (gamma) are average of the combustion chamber value and exit value
        return Thrust_mdot_cnst,Thrust,M_e

    def p_c_drop(self,chamber_ablation_rate, burnTime,showplot=True):
        erate = chamber_ablation_rate #m/s
        burnTimes = np.arange(0,burnTime,.01) #SECONDS
        p1 = self.engine.P_inj*1e5 # Pascals!
        ri = self.engine.engineProps[0, 0] # Meters!
        l = self.engine.engineProps[self.engine.throatInd,1] #Meters

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
                p2 = np.append(p2,p2_val)

            t = burnTimes
            name = ['Erate: ' + str(erate)]
            plt.figure()
            plt.show()
            plt.plot(t, p2/6894.76, label=name)
            plt.xlabel('Time (s)')
            plt.ylabel('Chamber Pressure (psi)')


        return pf


# ablative wall thickness sizing
    def size_ablator_thickness(self, burn_time,T_w, T_amb=272,wall_idx=[],t_A_design=[]):
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

        #Material properties:
        liner = self.SilicaPhenolic()
        OW = self.CarbonEpoxy()

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

    def size_throat_insert_thickness(self, burn_time,T_w, T_amb=272,wall_idx=[],t_C=.4):
        # Graphite thickness (m)
        t_A = t_C

        # Tamb is # 30F since colder is more conservative
        if wall_idx == []:
            wall_idx = self.engine.throatInd

        #Material properties:
        liner = self.Graphite()
        OW = self.CarbonEpoxy()

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




    def wall_temp_plot(self,T_w,mat,t,T_amb=272,numpts=200,wall_idx=[]):

        if wall_idx == []:
            wall_idx = self.engine.throatInd
        (h_g, q_conv, T_aw) = bartz(self.engine, T_w, wall_idx)


        k = mat.k


    def chamber_stress_analysis(self,material,t,DT):
        p_c = self.engine.P_inj*100000
        r_c = self.engine.engineProps[0, 0]
        CTE = material.CTE
        E = material.E
        v = material.v
        UTS=material.UTS

        stress_Hoop = p_c*r_c/t #Pa
        stress_Thermal =2*CTE*E*DT/(1-v)
        stress = stress_Hoop+stress_Thermal
        FOS = UTS/stress
        return stress,FOS


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
        Fo_arr = mat.k / mat.rho / mat.Cp * time_arr / thickness**2
        T_arr = np.zeros(len(x_star_arr))
        for i in np.arange(len(x_star_arr)):
            if i == 0:
                for z in np.arange(3):
                    Bi = h_g * thickness / mat.k
                    T = self.transient_wall_temp(Bi, Fo_arr, self.T_wall_i, self.engine.engineProps[station, 9], x_star_arr[i])[T]
                    h_g = bartz(self.engine, np.mean(T) , station)[0]
            else:
                Bi = h_g * thickness / mat.k
                T = self.transient_wall_temp(Bi, Fo_arr, self.T_wall_i, self.engine.engineProps[station, 9], x_star_arr[i])[T]
            T_arr[i] = T[-1]
        T_arr = np.flip(T_arr)

        plt.figure()
        plt.plot(x_star_arr * thickness / .0254, T_arr)
        plt.xlabel('Thickness (in)')
        plt.ylabel('Wall Temp (K)')
        plt.title('Temperature across Wall Thickness, Chamber')
        plt.show()





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


