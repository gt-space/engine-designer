import numpy as np
import scipy as sci
from matplotlib import pyplot as plt

import sys
# Uncomment and modify the following sys.path.insert line with the path of your file if "ModuleNotFoundError: No module named '...'" exception is given:
sys.path.insert(0, '/Users/jacob/Documents/Projects/Research/Film Cooling/engine-designer-master/engineDesigner_v5.0')

from regenDesigner.bartz import bartz
from regenDesigner.bartz import hg_gas_film
from regenDesigner.bartz import hg_boiling_liquid_film
from liquid_film_cooling import liquid_film_cooling
from fuel_props import JetA

class Heatsink:
    def __init__(self, engine, thickness, hotfire_time, chamber_inner_diameter, lfc, graphite_OD = 0.75, graphite_start_index = 130, graphite_end_index = 159, dt = 1.0, radial_subdivisions = 1000, T_initial = 298.15): #, axial_subdivisions = 10):
        # Eagle, H. A., & Ross, A. L. (1959). 
        # Steady-state and transient thermal stresses in a tube subjected to internal heating (Report No. APEX 460). 
        # General Electric Company, Atomic Products Division. 
        # https://doi.org/10.2172/4240259   
        self.engine = engine # Engine object from design script
        self.thickness = thickness # in
        self.Ri = chamber_inner_diameter/2 # in
        #self.graphite_thickness = graphite_thickness # in, radial thickness, assumed to be constant
        self.graphite_start_index = graphite_start_index # in, axial length of insert
        self.graphite_end_index = graphite_end_index # in, starting distance along axis from injector
        self.graphite_OD = graphite_OD # in, graphite insert outer diameter

        self.hotfire_time = hotfire_time # s
        self.dt = dt # s, transient time-step
        self.rad_sub = radial_subdivisions # radial divisions for finite difference method
        self.time_list = np.arange(0, hotfire_time+dt, dt) # list of all of the times that the FDM analysis is conducted at 
        self.thickness_array = (self.Ri - self.engine.engineProps[:, 0] * 39.3701) + self.thickness # in, indicates wall thickness at each axial station
        self.material_distribution = self.material()[0] # given the material assignment at each FDM analsyis point along the contour
        self.T_initial = T_initial # K
        self.temps = np.full((int(np.ceil(hotfire_time/dt)) + 1, len(self.engine.engineProps[:, 1]), self.rad_sub + 2), self.T_initial) # K, starting temperature distribution
        
        self.mat = self.material()[1]
        # self.axial_sub = axial_subdivisions # axial stations where thermal gradient is evaluated

        if len(engine.film_cooling) > 0:
            self.liquid_hg = self.lfc_hg(lfc, engine.film_cooling[5])

# -------- Material Definitions --------- #
# if a new material is defined, ensure that all of the properties are defined for each material:
# 1. thermal conductivity (k)
# 2. specific heat capacity (c)
# 3. density (rho)
# 4. melting temperature

    class Steel:
        def __init__(self):
            self.yield_strength = 370 # MPa, yield strength
            self.E = 205 # GPa, Young's Modulus
            self.c = 486 # J/kg-K, specific heat
            self.k = 51.9 # W/m-K, thermal conductivity
            self.rho = 7870 # kg/m^3, density
            self.CTE = 10e-6 # m/m-K, linear coefficient of thermal expansion
            self.v = 0.29 # Poisson's Ratio
            self.melting_temp = 1738 # K

    class Graphite:
        def __init__(self):
            #self.yield_strength = 31.9 # ksi, yield strength
            self.E = 27.6 # GPa, Young's Modulus
            self.CTE = 2e-6 # m/m-K, linear coefficient of thermal expansion
            self.v = 0.31 # Poisson's Ratio
            self.k = 140  # W/m-K
            self.rho = 1760  # kg/m^3
            self.c = 707 # J/kg-K

# -------- FDM Analysis -------- # 
    def material(self):
        # Generates the material distribution of steel and graphite given in self.material
        axial = self.engine.engineProps[:, 1] * 39.37 #  in
        material_distribution = np.full((len(self.engine.engineProps[:, 1]), self.rad_sub + 2), self.Steel())
        mat_string = np.full((len(self.engine.engineProps[:, 1]), self.rad_sub + 2), 'S')
        for z in range(len(axial)):
            if z >= self.graphite_start_index and z <= self.graphite_end_index:
                g_thick = (self.graphite_OD/2) - (self.engine.engineProps[z, 0] * 39.3701)
                p = np.floor((g_thick/self.thickness_array[z]) * (self.rad_sub + 1))
                for r in range(int(p)):
                    material_distribution[z, r] = self.Graphite()
                    mat_string[z, r] = 'g'
        return material_distribution, mat_string 
    # material_distribution contains array of material objects
    # mat_string contains string array of materials as printable string

    def transient(self, hg, t, z):
        # Calculates the temperature distribution at a given time in time_list and axial station by using transient FDM matrix given in https://doi.org/10.2172/4240259
        init_temps = np.array(self.temps[t-1, z, :])
        a = self.engine.engineProps[z, 0] # m, inner radius
        b = (self.thickness + self.Ri) / 39.37 # m, outer radius
        l = np.log(b/a) # non-dimensionalized thickness
        d_eta = l/self.rad_sub # non-dimensionalized differential in radial subdivisions
        eta = np.linspace(0, l + d_eta, self.rad_sub + 2)
        inner_mat = self.material_distribution[z, 0]
        k = inner_mat.k # W/m-K
        phi = k / (a * hg * d_eta)
        k_array = np.array([i.k for i in self.material_distribution[z, :]]) # W/m-K
        rho_array = np.array([i.rho for i in self.material_distribution[z, :]]) # kg/m^3
        c_array = np.array([i.c for i in self.material_distribution[z, :]]) # J/kg-K
        Ko = (c_array * rho_array * a**2 * np.exp(2 * eta) * d_eta**2) / (k_array * self.dt)
        #Ko = c_array
        Co = -1*Ko*init_temps
        Co[-1] = 0
        Co[0] = self.engine.engineProps[z, 9] # K
        x = -1*(2+Ko)
        x[0] = 1 + phi
        x[-1] = -1
        Ao = np.diag(x)
        for i in range(0, len(x)-1):
            Ao[i, i-1] = 1
            Ao[i, i+1] = 1
        Ao[0, 1] = -1*phi
        Ao[0, -1] = 5
        Ao[-1, -3] = 1
        final_temps = np.matmul((np.linalg.inv(Ao)), Co)
        self.temps[t, z, :] = final_temps
        return init_temps, final_temps
    
    def iterate(self,t,z, film_cool_type, last_wall_temp, last_u_cool=None,gas_film_present=None):
        # Finds wall temperature and convective coeffcient by iterating Bartz correlation several times at each axial station for a given time in time_list
        dz = (self.engine.engineProps[1,1]-self.engine.engineProps[0,1])
        T_hg = self.engine.engineProps[z,9] # approximate mainstream gas temp to be gas temp without film cooling
        u_cool = -np.inf
        heat_flux_wall = .0 # can update later to get more accurqate h_g for liquid film
        deltaQ = 0 # can update later to get more accurate u_cool values for gas film
        #print(f'z:{z},gas film present?: {gas_film_present}')
        #print(self.engine.film_cooling[-1][2])
        for i in range(0, 10):
            if film_cool_type == "liquid" and (z * dz) < self.engine.film_cooling[-1][2] and z > 5:
                M_wt = self.engine.film_cooling[-1][1]
                #hg  = hg_boiling_liquid_film(self.engine, self.temps[t, z, 0], T_hg, heat_flux_wall, M_wt, z, dz, lfc)
                hg = self.liquid_hg
            elif film_cool_type == "gas" and gas_film_present and z > 5:
                M_wt = self.engine.film_cooling[-1][1]
                hg, u_cool, gas_film_present = hg_gas_film(self.engine, self.temps[t, z, 0],T_hg, deltaQ, last_wall_temp, last_u_cool, M_wt, z, dz)
            else:
                hg = bartz(self.engine, self.temps[t, z, 0], z)[0]
            _, final_temps = self.transient(hg, t, z)
            self.temps[t, z, :] = final_temps
        return hg, final_temps[0], u_cool, gas_film_present
    
    @staticmethod
    def lfc_hg(param, axial_velocity_fraction):
        _, lfc_hg, _ = liquid_film_cooling.get_film_cooled_length(param, axial_velocity_fraction)
        return lfc_hg

    
    # def transient_solution(self):
    #     # Iterates along every axial station and time in time_list to find overall temperature history at each axial station along the wall thickness
    #     hg_list = []
    #     last_wall_temp = 298.15 # (K)
    #     if len(self.engine.film_cooling) > 0:
    #         last_u_cool = self.engine.film_cooling[-1][-1]
    #         film_cool_type = self.engine.film_cooling[0]
    #     else:
    #         last_u_cool = -np.inf
    #         film_cool_type = None
    #         gas_film_present = True # if engine is film cooled with gas, the coolant is separate from mainstream flow at the injector
    #     for z in range(len(self.temps[0, :, 0])):
    #         for t in range(1, len(self.temps[:, 0, 0])):
    #             hg, last_wall_temp, last_u_cool, gas_film_present = self.iterate(t,z,film_cool_type=film_cool_type,last_wall_temp=last_wall_temp,last_u_cool=last_u_cool, gas_film_present=gas_film_present)
    #         hg_list.append(hg)
    #     if len(self.engine.film_cooling) > 0: # and not gas_film_present:
    #         self.engine = self.engine.film_cooling[-2] # replace Engine with updated MR Engine if homogeneous temp. is reached
    #     return hg, final_temps

    def transient_solution(self):
        # Iterates along every axial station and time in time_list to find overall temperature history at each axial station along the wall thickness
        hg_list = np.zeros((len(self.time_list), len(self.engine.engineProps[:, 1])))
        for z in range(len(self.engine.engineProps[:, 1])):
            for t in range(1, len(self.time_list)):
                hg, _, _, _ = self.iterate(t, z,  self.engine.film_cooling[0] if len(self.engine.film_cooling) > 0 else "", last_wall_temp=298.15)
                hg_list[t, z] = hg
        return self.temps, hg_list
    
    def closest(self, time):
        # Finds the closest time in time_list to the given time

        #length_list = self.engine.engineProps[:, 1] * 39.37 # m --> in
        #z_index = np.argmin([np.abs(l - length) for l in length_list]) # gets closest length value in temperature list in inches
        t_index = np.argmin([np.abs(x - time) for x in self.time_list]) # gets closest time value in temperature list in seconds 
        return t_index #, z_index 
    
    def r_array(self, z):
        # Returns an array of radial indices at a given axial station
        #_, z = self.closest(0, length)
        a = self.engine.engineProps[z, 0] # m, inner radius
        b = (self.thickness + self.Ri) / 39.37 # m, outer radius
        l = np.log(b/a) # non-dimensionalized thickness
        d_eta = l/self.rad_sub # non-dimensionalized differential in radial subdivisions
        eta = np.linspace(0, l + d_eta, self.rad_sub + 2)
        r = (np.exp(eta) * a * 39.37) - (self.engine.engineProps[z, 0] * 39.37) # in, radial locations
        return a, b, r

    def wall_temp_2d(self, time, z):
        # Finds wall temperature along thickness at a given axial station and time
        _, hg_list = self.transient_solution()
        t = self.closest(time)
        a, b, r = self.r_array(z)
        p = np.polyfit(r, self.temps[t, z, :], 5)
        r = np.linspace(0, self.thickness, 10)
        T = np.polyval(p, r)
        hg = hg_list[z]
        return a, b, hg, r, T,

    def plot_wall_temp_gradient_at_station(self, time, z):
        # Plots wall temperature gradient along thickness at a given axial station and time
        _, _, _, r, T = self.wall_temp_2d(time, z)
        x = self.engine.engineProps[z, 1] * 39.37
        plt.figure()
        plt.plot(r, T)
        plt.xlabel('Thickness (in)')
        plt.ylabel('Wall Temp (K)')
        plt.title(f'Temperature across Wall Thickness @ {x:.1f} in. and {time: .1f} sec.\n Max Temperature: {T[0]:.1f} K')
        return r, T
    
    def plot_inner_wall_temp_at_time(self, time): 
        # Plots inner wall temperature at a given time along entire contour
        self.transient_solution()
        t = self.closest(time)
        plt.figure()
        plt.plot(self.engine.engineProps[:, 1] * 39.37, self.temps[t, :, 0])
        plt.xlabel('Axial Location (in)')
        plt.ylabel('Inner Wall Temp (K)')
        plt.title(f'Temperature along Inner Wall @ {self.time_list[t]:.1f} sec.')
    
    def plot_transient_3d(self):
        # Plots transient temperature distribution along entire contour
        self.transient_solution()
        t,Z = np.meshgrid(self.time_list[1:-1], self.engine.engineProps[:, 1] * 39.37)
        plt.figure()
        ax_ch = plt.axes(projection='3d')
        ax_ch.plot_surface(Z, t, np.transpose(self.temps[1:-1, :, 0]), rstride=1, cstride=1,
                               cmap='viridis', edgecolor='none')
        ax_ch.set_xlabel('Axial Location (in)')
        ax_ch.set_ylabel('Time (sec)')
        ax_ch.set_zlabel('Wall Temp (K)')
        plt.title('Wall Temperature over time')
        #return t, self.temps[:, :, 0]

    def plot_inner_wall_temp_at_station(self, z):
        # Plots inner wall temperature at a given axial station over time
        self.transient_solution()
        x = self.engine.engineProps[z, 1] * 39.37
        plt.figure()
        plt.plot(self.time_list, self.temps[:, z, 0])
        plt.xlabel('Time (sec)')
        plt.ylabel('Inner Wall Temp (K)')
        plt.title(f'Temperature at Inner Wall @ {x:.2f} in. over time')

    def graphite_test(self, temp_variation=100):
        # Uses melting temp in given steel properties to find starting and ending axial indices of graphite insert. Compares to current chosen endpoints for insert
        test_mat = np.full((len(self.engine.engineProps[:, 1]), self.rad_sub + 2), self.Steel())
        test_temps = np.full((int(np.ceil(self.hotfire_time/self.dt)) + 1, len(self.engine.engineProps[:, 1]), self.rad_sub + 2), self.T_initial)

        def test_transient(self, hg, t, z):
            init_temps = np.array(test_temps[t-1, z, :])
            a = self.engine.engineProps[z, 0] # m, inner radius
            b = (self.thickness + self.Ri) / 39.37 # m, outer radius
            l = np.log(b/a) # non-dimensionalized thickness
            d_eta = l/self.rad_sub # non-dimensionalized differential in radial subdivisions
            eta = np.linspace(0, l + d_eta, self.rad_sub + 2)
            inner_mat = test_mat[z, 0]
            k = inner_mat.k # W/m-K
            phi = k / (a * hg * d_eta)
            k_array = np.array([i.k for i in test_mat[z, :]]) # W/m-K
            rho_array = np.array([i.rho for i in test_mat[z, :]]) # kg/m^3
            c_array = np.array([i.c for i in test_mat[z, :]]) # J/kg-K
            Ko = (c_array * rho_array * a**2 * np.exp(2 * eta) * d_eta**2) / (k_array * self.dt)
            #Ko = c_array
            Co = -1*Ko*init_temps
            Co[-1] = 0
            Co[0] = self.engine.engineProps[z, 9] # K
            x = -1*(2+Ko)
            x[0] = 1 + phi
            x[-1] = -1
            Ao = np.diag(x)
            for i in range(0, len(x)-1):
                Ao[i, i-1] = 1
                Ao[i, i+1] = 1
            Ao[0, 1] = -1*phi
            Ao[0, -1] = 0
            Ao[-1, -3] = 1
            final_temps = np.matmul((np.linalg.inv(Ao)), Co)
            test_temps[t, z, :] = final_temps
            return init_temps, final_temps
        
        def test_iterate(self, t, z):
            final_temps = test_temps[t, z, :]
            for i in range(0, 10):
                hg = bartz(self.engine, test_temps[t, z, 0], z)[0]
                _, final_temps = test_transient(self, hg, t, z)
                test_temps[t, z, :] = final_temps
            return hg, final_temps
        
        def test_transient_solution(self):
            for z in range(len(test_temps[0, :, 0])):
                for t in range(1, len(self.temps[:, 0, 0])):
                    test_iterate(self, t, z)
            return test_temps
    
        test_temps = test_transient_solution(self)[-1, :, 0]

        graphite_list = []

        for z in range(len(test_temps)):
            if test_temps[z] >= self.Steel().melting_temp - temp_variation:
                graphite_list.append(z)

        test_graphite_start = graphite_list[0]
        test_graphite_end = graphite_list[-1]
        end_min = int(np.floor(len(self.engine.engineProps[:, 0])*2/3))
        z_min = np.argmin([np.abs(self.graphite_OD/2 - r) for r in (self.engine.engineProps[0:end_min, 0] * 39.37)])
        min_thick = self.graphite_OD/2 - self.engine.engineProps[self.graphite_start_index, 0]*39.37
        return f'Start graphite at index {test_graphite_start} ({self.engine.engineProps[test_graphite_start, 1]*39.37:.2f} in.) and end at index {test_graphite_end} ({self.engine.engineProps[test_graphite_end, 1]*39.37:.2f} in.).\nMinimum starting index is {z_min}.\nCurrently starting at index {self.graphite_start_index} ({self.engine.engineProps[self.graphite_start_index, 1]*39.37:.2f} in.) and ending at {self.graphite_end_index} ({self.engine.engineProps[self.graphite_end_index, 1]*39.37:.2f} in.)\nCurrent thickness: {min_thick:.2f} in\n '

# Thermal Stress values are not accurate: will be fixed in future update
'''
    def thermal_stresses(self, time, z): # at inner wall (where stresses are max)
        a, b, r, T = self.wall_temp_2d(time, z)
        a = a * 39.37
        b = b * 39.37
        dr = r[1] - r[0]
        t = self.closest(time)
        r = r + (self.engine.engineProps[z, 0] * 39.37)
        T_av = (2 / (b**2 - a**2)) * np.sum(T*r*dr) # K, relative avg. temp to ambient
        T_inner = T[0] - self.T_initial # K, relative inner wall temp to ambient
        inner_mat = self.material_distribution[z, 0]
        cte = inner_mat.CTE # m/m-K
        E = inner_mat.E # GPa
        v = inner_mat.v
        von_mises = np.abs((cte*E/(1-v)) * (T_av - T_inner)) # hoop is same axial stress at inner wall, radial thermal stress is zero
        return T_av'''


'''   def transient_solution(self):
        # Iterates along every axial station and time in time_list to find overall temperature history at each axial station along the wall thickness
        hg_list = []
        for z in range(len(self.temps[0, :, 0])):
            for t in range(1, len(self.temps[:, 0, 0])):
                hg, _ = self.iterate(t, z)
            hg_list.append(hg)
        return self.temps, hg_list
    # temps: full wall temperature history (3d array)
    # hg_list: list of convective coefficients at all axial stations'''