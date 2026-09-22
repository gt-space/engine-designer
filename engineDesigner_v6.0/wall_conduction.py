def wall_temps(film, bulk_velo, x, Taw, copper_k, wall_thickness,backside_temp, backside_Q):

    x_safe = max(x, 1 * (10**-5))

    #local reynolds number
    reynolds = (film.film_rho * bulk_velo * x_safe) / film.mu

    #chilton-colburn correlation heat transfer equation
    nusselt = 0.0296 * (reynolds **0.8) * (film.film_Pr ** 0.3)
    #film thermal resistance
    h_film = (nusselt * film.k) / x_safe
    film_tr = 1 / h_film

    cu_tr = wall_thickness / copper_k
    backside_tr = 1 / backside_Q

    total_tr = film_tr + cu_tr + backside_tr
    total_td = Taw - backside_temp

    film2wall_q = total_td / total_tr

    hot_wall = Taw - (film2wall_q * film_tr)
    cold_wall = hot_wall - (film2wall_q * cu_tr)

    return film2wall_q, hot_wall, cold_wall

