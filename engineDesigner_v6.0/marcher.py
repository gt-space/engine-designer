import yaml
import math
import matplotlib.pyplot as plt
import thermo_propertiesV2 as thermo
import boundary_layer as bl
import heat_mass_transferV2 as hmt
import wall_conduction as wc
import contour_kinematics as ck

def physics(x, linear_mdot, geometry, h_current, config):

    r = geometry.local_radius(x)
    local_area =  math.pi * (r ** 2)
    local_circ = 2 * math.pi * r

    mdot_circum = linear_mdot / local_circ
    Pstag = config['engine']['pressure_chamber_bar'] * 100000
    Tstag = config['engine']['temperature_core_gas']
    phi = config['engine']['equivalence_ratio']
    core_stag = thermo.CoreGas(Pstag, Tstag, 0, phi)
    gamma = core_stag.gamma
    R = core_stag.R

    Plocal, Tlocal, local_velo, mach_local = ck.gas_state(local_area, geometry.throat_area, Pstag, Tstag, gamma, R)

    #cantera states
    T_inject = config['engine']['temperature_injection']
    scFilm = thermo.SupercriticalFluid(Plocal, T_inject, H_current=h_current)
    core_gas = thermo.CoreGas(Plocal, Tlocal, local_velo, phi)

    #shear stress
    skin_fric_o, St_o, g_m, shear = bl.skin_friction(core_gas, x, r)
    
    #couette flow assumption for film velocity
    if linear_mdot <= 1 * (10 ** -8):
        return None #depleted film

    film_thickness = math.sqrt((2 * scFilm.mu * mdot_circum) / (scFilm.film_rho * shear))
    film_velo = mdot_circum / (scFilm.film_rho * film_thickness)

    #heat flux
    q_gas2film, diffusion_G, Taw = hmt.calc_mixing_and_gas_heat(core_gas, scFilm, St_o, g_m, mdot_circum, x)

    #wall conduction
    q_film2wall, hot_wall_temp, cold_wall_temp = wc.wall_temps(scFilm, film_velo, x, Taw, 
                                                               config['material']['thermal_conductivity_copper'], 
                                                               config['material']['thickness_copper_wall'], 
                                                               config['material']['temperature_backside_coolant'], 
                                                               config['material']['heat_transfer_coefficient_backside'])

    #differentials
    derivative_m = -(diffusion_G * local_circ)

    q_net2film = q_gas2film - q_film2wall
    derivative_h = q_net2film / mdot_circum

    return {'derivative_m': derivative_m, 'derivative_h': derivative_h, 'x': x, 'film_temp': scFilm.film_temp, 
            'film_H': scFilm.film_H,  'hot_wall_temp': hot_wall_temp, 'cold_wall_temp': cold_wall_temp, 'Taw': Taw, 
            'mdot_circum': mdot_circum, 'diffusion_G': diffusion_G, 'q_gas2film': q_gas2film, 'q_film2wall': q_film2wall, 
            'film_velo': film_velo, 'mach_local': mach_local}

def rk4_marcher(config_file, contour_file):

    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    geometry = ck.Contour(contour_file)

    step = config['solver']['step_size_axial']
    chamber_length = config['solver']['length_chamber_total']
    wall_temp_allow = config['solver']['max_allowable_wall_temp']

    x = 0
    P_initial = config['engine']['pressure_chamber_bar'] * 100000
    T_inject = config['engine']['temperature_injection']
    r = geometry.local_radius(0)
    linear_mdot = config['engine']['mass_flow_per_circumference'] * 2 * math.pi * r

    film_state1 = thermo.SupercriticalFluid(P_initial, T_initial=T_inject, T_inject=T_inject)
    h_current = film_state1.film_H

    data = []

    while x < chamber_length:

        k1 = physics(x, linear_mdot, geometry, h_current, config)
        if k1 is None:
            print(f"Film depleted at x = {x:.4f} m")
            break
        if k1['hot_wall_temp'] > wall_temp_allow:
            print(f"Wall temperature exceeded at x = {x:.4f} m")
            break

        data.append(k1)

        mid_x = x + step / 2
        mid_m_k1 = linear_mdot + (k1['derivative_m'] * (step / 2))
        mid_h_k1 = h_current + (k1['derivative_h'] * (step / 2))

        k2 = physics(mid_x, mid_m_k1, geometry, mid_h_k1, config)

        mid_m_k2 = linear_mdot + (k2['derivative_m'] * (step / 2))
        mid_h_k2 = h_current + (k2['derivative_h'] * (step / 2))

        k3 = physics(mid_x, mid_m_k2, geometry, mid_h_k2, config)

        end_x = x + step
        end_m_k3 = linear_mdot + (k3['derivative_m'] * step)
        end_h_k3 = h_current + (k3['derivative_h'] * step)

        k4 = physics(end_x, end_m_k3, geometry, end_h_k3, config)

        linear_mdot += (step / 6) * (k1['derivative_m'] + 2 * k2['derivative_m'] + 2 * k3['derivative_m'] + k4['derivative_m'])
        h_current += (step / 6) * (k1['derivative_h'] + 2 * k2['derivative_h'] + 2 * k3['derivative_h'] + k4['derivative_h'])

        x += step

    return data

def plotter(data):
    x_vals = [entry['x'] for entry in data]
    film_temp_vals = [entry['film_temp'] for entry in data]
    hot_wall_temp_vals = [entry['hot_wall_temp'] for entry in data]
    cold_wall_temp_vals = [entry['cold_wall_temp'] for entry in data]
    Taw_vals = [entry['Taw'] for entry in data]
    mdot_circum_vals = [entry['mdot_circum'] for entry in data]
    diffusion_G_vals = [entry['diffusion_G'] for entry in data]
    q_gas2film_vals = [entry['q_gas2film'] for entry in data]
    q_film2wall_vals = [entry['q_film2wall'] for entry in data]
    mach_local_vals = [entry['mach_local'] for entry in data]

    fig, axs = plt.subplots(3, 2, figsize=(14, 12))
    fig.suptitle('Film Cooling Simulation Results', fontsize=16)

    #temperature plot
    axs[0, 0].plot(x_vals, Taw_vals, label='Adiabatic Wall Temperature (K)', color='blue')
    axs[0, 0].plot(x_vals, hot_wall_temp_vals, label='Hot Wall Temperature (K)', color='red')
    axs[0, 0].plot(x_vals, cold_wall_temp_vals, label='Cold Wall Temperature (K)', color='green')
    axs[0, 0].plot(x_vals, film_temp_vals, label='Film Bulk Temperature (K)', color='orange')
    axs[0, 0].set_title('Temperature Profiles')
    axs[0, 0].set_xlabel('Axial Position (m)')
    axs[0, 0].set_ylabel('Temperature (K)')
    axs[0, 0].legend()
    axs[0, 0].grid()

    #coolant mass
    axs[2, 1].plot(x_vals, mdot_circum_vals, label='Coolant Mass Flow Rate (kg/s/m)', color='purple')
    axs[2, 1].set_title('Coolant Mass Flow Rate')
    axs[2, 1].set_xlabel('Axial Position (m)')
    axs[2, 1].set_ylabel('Gamma (kg/ m*s)')
    axs[2, 1].grid()

    #heat fluxes
    axs[1, 0].plot(x_vals, q_gas2film_vals, label='Heat Flux from Gas to Film (W/m^2)', color='cyan')
    axs[1, 0].plot(x_vals, q_film2wall_vals, label='Heat Flux from Film to Wall (W/m^2)', color='magenta')
    axs[1, 0].set_title('Heat Fluxes')
    axs[1, 0].set_xlabel('Axial Position (m)')
    axs[1, 0].set_ylabel('Heat Flux (MW/m^2)')
    axs[1, 0].legend()
    axs[1, 0].grid()

    #mass diffusion
    axs[1, 1].plot(x_vals, diffusion_G_vals, label='Mass Diffusion Rate (kg/s/m^2)', color='brown')
    axs[1, 1].set_title('Mass Diffusion Rate')
    axs[1, 1].set_xlabel('Axial Position (m)')
    axs[1, 1].set_ylabel('Mass Lost to Core Gas (kg/ m^2 *s)')
    axs[1, 1].grid()

    #core gas
    axs[2, 0].plot(x_vals, mach_local_vals, label='Core Gas Mach Number', color='black')
    axs[2, 0].set_title('Core Gas Mach Number')
    axs[2, 0].set_xlabel('Axial Position (m)')
    axs[2, 0].set_ylabel('Mach Number')
    axs[2, 0].grid()

    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    plt.show()

if __name__ == "__main__":
    config_file = 'config1.yaml'
    contour_file = 'contour1.csv'

    data = rk4_marcher(config_file, contour_file)
    plotter(data)





