import math

def calc_mixing_and_gas_heat(gas, film, St_o, g_m, mdot_circum, x):
    x_safe = max(x, 1 * (10**-5))

    #spalding mass transfer number = force that pulls coolant into core gas
    cp_avg = (gas.cp + film.film_cp) / 2
    spalding = (cp_avg * (gas.temp - film.film_temp)) / film.q_barrier

    #mass diffusion rate
    diffusion_mdot = g_m * math.log(1+spalding)

    #transpiration blowing reduction = mass that mixes outward shields wall
    if spalding > 0:
        correction = math.log(1+spalding) / spalding
    else:
        correction = 1.0

    St_blow = St_o * correction

    #heat flux from core gas going into film
    gas2film_H = gas.core_velo * gas.rho * gas.cp * St_blow
    gas2film_Q = gas2film_H * (gas.temp - film.film_temp)

    #hatch papell for adiabatic wall temp
    film_TC = max(mdot_circum * film.film_cp, 1)
    mixing_exp = math.exp(-(gas2film_H * x_safe) / film_TC)
    Taw = gas.temp - (gas.temp - film.film_temp) * mixing_exp

    return gas2film_Q, diffusion_mdot, Taw



