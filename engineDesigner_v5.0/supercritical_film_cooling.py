# assume that coolant becomes supercritical upon injection

import numpy as np

class supercritical_film_cooling:

    def __init__(self, cea_obj, mdot_gas0, MR, contour, eps,
                 mdot_cool, pressure_orifice_cool, temp_orifice_cool, d_film_orifice, num_film_orifices, orifice_cstar):
        
        self.cea_obj = cea_obj
        self.mdot_gas0 = mdot_gas0
        self.MR = MR
        self.contour = contour
        self.eps = eps
        self.mdot_cool = mdot_cool
        self.pressure_orifice_cool = pressure_orifice_cool
        self.temp_orifice_cool = temp_orifice_cool
        self.d_film_orifice = d_film_orifice
        self.num_film_orifices = num_film_orifices
        self.orifice_cstar = orifice_cstar
    
    def solve(self, num_nodes, tol, max_counts):
        node_volume = self.engine_volume / num_nodes
        delta_x = node_volume**(1/3)
        # guess h_old ; probably use cea or other part of the engine code
        counter = 0
        while True:
            (h, wall_temps, film_cooled_length) = self.iterate(h_old, num_nodes, node_volume, delta_x)
            if abs(h-h_old)<tol: # check for convergence 
                return (h, wall_temps, film_cooled_length)
            h_old = h
            counter += 1
            if counter > max_counts:
                print("Film cooling solver failed to converge")
                return (h, wall_temps, film_cooled_length)
    
    # nodes are cubes 
    # assume values for density, temperature, axial velocity (u), and radial velocity (v)
    # effing no
    def iterate(self, h, num_nodes, node_volume, delta_x):
        wall_temps = []
        dz = 1e-3 # meters


        # use gas governing equations to obtain internal wall temp ; add this to temp list

        # A is a 6 x 6 matrix with a row for each governing equation (each column for a different node)
        # X is a 6 x 6 matrix with a row for each node, each column contains the a properties of the node centroid
        # M is a 6 x 6 matrix, with each row being the RHS of each governing equation, and each RHS a different node
        # properties: rho (kg/m^3), u (axial speed, m/s), v (radial speed, m/s), pressure (Pa), temp (K), concentration of film
        A = [[]]
        M = [[]]
        X = [[]]

        rho_matrix = [[]]

        # loop through all of the nodes
        z = 0 # axial position (m), start at injection
        for i in node_groups:
            r = 0 # radial position (m), start at center
            for j in np.linspace(1, node_groups(i)):

                # first governing equation: conservation of mass --> solve for rho
                nextr_rho_coeff_RHS = -last_r_mass_conv_flux - last_z_mass_conv_flux # flux from other nodes
                nextr_rho_coeff_RHS -= area * 
                nextr_rho_coeff = area*
                lastr_rho_coeff = 1/r*lastr_r*lastr_v(del_r)
                nextz_rho_coeff = nextz_u/del_z
                lastz_rho_coeff = lastz_u/del_z
                rho_RHS = np.zeros(2, 1)
                rho_matrix[0] = [nextr_rho_coeff, lastr_rho_coeff]
                rho_matrix[1] = [nextz_rho_coeff, lastz_rho_coeff]
                rho_matrix = 

                # second governing equation: ideal gas law --> solve for pressure
                pressure = rho*R*temp


                # third governing equation: conservation of axial momentum --> solve for u
                nextz_u_coeff = rho*u/delta_x - eta

                # fourth governing equation: conservation of radial momentum --> solve for v
                nextz_temp_coeff = 

                # fifth governing equation: concentration equation --> solve for c

                # fifth governing equation: conservation of energy --> solve for temperature
                
                
                if i is not node_groups(i):
                    # update to next radial position
                    r += delta_x
                    rho = lastr_rho
                    u = lastr_u
                    v = lastr_v
                    p = lastr_p
                    temp = lastr_temp
                    c = lastr_c
                else:
                    h.append(q_out/(area*(temp-self.temp_outside)))
            z += delta_x # update axial position

       # this is

