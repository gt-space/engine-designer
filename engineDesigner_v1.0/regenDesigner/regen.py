import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from .fuel_props import get_props
from .bartz import bartz
# from .pressureDrop import pressureDrop

np.set_printoptions(threshold=sys.maxsize) # Print full arrays (for debugging)

class RegenJacket:
    # Initialize the jacket object upon declaration
    def __init__(self, engine, channel_h=0.002, wall_t=0.002, min_fin_w = 0.001, channel_w = 0.002, T_ci=300, start_ind = 40):
        self.engine = engine # Engine object to be jacketed
        self.channel_h = channel_h # Channel height [m]
        self.wall_t = wall_t # Inner wall thickness [m]
        self.min_fin_w = min_fin_w # Minimum desired fin width [m]
        self.channel_w = channel_w # Width of channel [m]
        self.T_ci = T_ci # Coolant inlet temperature [K]
        self.start_ind = start_ind # Index where fuel enters (0 is exit plane)

    # Friction correlation function
    def get_friction(self, Re, D_hyd):
        def colebrook(Re, D_hyd):
            import numpy
            from scipy.optimize import root

            k = 0.005 * (10 ** (-3)) #Surface roughness, guess based on our manufacuting capability, mm => m

            def f(x):
                return -2*numpy.log10(k/(3.71*D_hyd) + 2.51/(Re*numpy.sqrt(x))) - 1.0/numpy.sqrt(x)

            fric = root(f, 0.04).x[0]
            return fric

        # Solve for Darcy Friction Factor Approximation using Reynold's number
        if Re <= 2320:
            fric = 64 / Re # Laminar flow approximation
        else:
            fric = colebrook(Re, D_hyd) # Colebrook-White Equation solution
        return fric

    # Simulate results given the geometry
    def simulate_regen(self):

        # ==== Define properties ====
        mDot = self.engine.mDot_f # Mass flow of coolant [kg/s]
        cond_w = 330 # copper conductivity [W/m-K]
        T_cb = self.T_ci # Initial bulk temp [K]
        P_co = self.engine.P_inj * 1.2 # Initial coolant pressure [bar] 1.2 included as stiffness

        # Initialize output vectors
        numPTS = np.size(self.engine.engineProps,0) # Number of engine stations
        wall_temps = np.zeros([numPTS-1-self.start_ind,1]) # Wall temperature vector
        coolant_temps = np.zeros([numPTS-1-self.start_ind,1]) # Coolant temperature vector
        pressures = np.zeros([numPTS-1-self.start_ind,1]) # Coolant pressure vector

        num_channels = round(2 * math.pi * (self.engine.R_t + self.wall_t) / (self.min_fin_w + self.channel_w))
        print(num_channels)

        # Main loop
        for x in range(numPTS-1-self.start_ind):
            i = numPTS - x - self.start_ind - 1 # Reverse order
            T_wg = 800 # Starting gas side wall temp estimate
            q_in = 1 # bs initlialization
            q_out = 2 # bs initlialization
            # Initialize parameters
            (rho_c, C_pc, cond_c, viscK_c) = get_props(T_cb) # Get coolant fluid properties
            R = self.engine.engineProps[i,0]
            fin_w = ((2 * math.pi * (R + self.wall_t)) - (num_channels * self.channel_w)) / num_channels
            length = self.engine.engineProps[i,1] - self.engine.engineProps[i-1,1] # Station Length [m]
            mDot_chan = mDot/num_channels

            while abs((q_in-q_out)/q_in) > 0.001:
                (h_g, qdot_ge, T_aw) = bartz(self.engine, T_wg, i)

                T_wc = T_wg - qdot_ge * self.wall_t/cond_w # Solve for coolant side wall temp

                D_hyd = 2 * self.channel_w * self.channel_h / (self.channel_w + self.channel_h) # Channel hydraulic diameter
                A_cc = self.channel_w * self.channel_h # Cross sectional area of channel (m^2)
                vel_c = mDot_chan/(A_cc * rho_c) # Coolant velocity in channel (m/s)
                Re_c = (vel_c * D_hyd) / viscK_c # Coolant Reynolds number
                Pr_c = viscK_c * rho_c * C_pc / cond_c # Coolant Prandtl Number

                # Apply Nusselt relation (see H&H pg 90):
                # h_c = 0.021 * (Re_c)**0.8 * (Pr_c)**0.4 * (0.64 + 0.36 * (T_cb/T_wc)) * (cond_c/D_hyd) # Old code
                viscK_c_w = get_props(T_wc)[3] # Get viscosity at wall
                # h_c = 0.0214 * Re_c**0.8 * (Pr_c)**0.4 * (1000000*viscK_c/(0.2))**0.14 * (cond_c/D_hyd)
                # h_c_2 = 0.021 * (Re_c**0.8) * (Pr_c**0.4) * (0.64+0.36*(T_cb/T_wc)) * (cond_c/D_hyd)
                h_c = 0.023 * (Re_c**0.8) * (Pr_c**0.4) * (cond_c/D_hyd) # Dittus Boelter reccomended for small dT Convection

                # Friction correlation (set up later)
                f = self.get_friction(Re_c, D_hyd)
                # f = (0.79 * np.log(Re_c) - 1.64) ** (-2)
                Nu = (f/8)*(Re_c - 1000)*Pr_c/(1 + 12.7*((f/8)**0.5)*(Pr_c**(2/3)-1))
                h_c = Nu * cond_c / D_hyd * (1000000*viscK_c/(0.2))**0.11


                # Future versions should edit this out of calculation
                #Run fin efficiency calculations to find the effective contact area
                # m = math.sqrt((h_c * perim_conf)/(cond_w * A_conf)) # Old fin efficiency parameter
                m = math.sqrt((2*h_c)/(cond_w*fin_w)) # Fin efficiency length independent
                eta_fin = math.tanh(m * self.channel_h)/(m * self.channel_h) # Fin efficiency
                A_fin = 2 * self.channel_h * length # Area of single fin in analysis segment
                A_finTot = num_channels * A_fin # Total fin area
                A_wallTot = num_channels * self.channel_w * length  # Total wall contact area
                A_totc = A_finTot + A_wallTot # Total coolant side heat transfer area
                eta_tot = (A_wallTot + A_finTot * eta_fin)/A_totc  # Overall convection efficiency
                q_out = h_c * A_totc * eta_tot * (T_wc - T_cb)

                A_totg = 2 * math.pi * R * length #total hot gas side wall area
                q_in = qdot_ge * A_totg # Total heat entering wall from gasses (not the rate)

                # If too much heat is entering
                if q_in > q_out:
                    T_wg = T_wg * 1.001
                else:
                    T_wg = T_wg * .999

            dT_c = ((q_in + q_out)/2)/(C_pc * mDot) # Temperature change in coolant
            dP = f * ((length)/D_hyd) * rho_c * ((vel_c**2)/2) / 100000 #Solve for final dP across passage segment (Pa > Bar)

            T_cb += dT_c
            P_co -= dP

            wall_temps[i-1,0] = T_wg
            coolant_temps[i-1,0] = T_cb
            pressures[i-1,0] = P_co

        plt.plot(self.engine.engineProps[:-(self.start_ind+1):,1], wall_temps)
        plt.xlabel('Distance from Injector (m)', fontsize=16)
        plt.ylabel('Coolant Temperature (K)', fontsize=16)
        plt.show()

        return (wall_temps, coolant_temps, pressures, num_channels)
