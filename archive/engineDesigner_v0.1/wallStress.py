# Script for calculating wall thickness given the expected stresses

def wallStress(P_ce_req, D_hyd, E, a, v, qdot_ge, cond_w, yield_str):
    # Run convergence to find wall_t
    wall_t = 0.001 #Initial thickness, m
    thicken = True
    stress_old = 0
    while True:
        # Run channel wall stress analysis
        stress_dP = (P_ce_req * 5 - 101325) * D_hyd/wall_t #5 is arbitrary FOS for now
        stress_Therm = E * a * qdot_ge * wall_t/(2 * (1 - v)*cond_w)
        stress_tot = stress_dP + stress_Therm
        # If the stress went up from thickening, flip to thinning and visa versa
        if stress_old < stress_tot:
            thicken = not thicken

        # For now to be safe we will set yielding as the limit
        if stress_tot > yield_str:
            # Either we're too thin (high P stress) or too thick (high T stress)
            # thicken the wall and see what happens on next iteration
            stress_old = stress_tot
            if thicken:
                wall_t *= 1.001
            else:
                wall_t *= 0.999
        # If we're below yield, shrink thickness
        else:
            wall_t *= 0.999
        # Check if we're at the stress and thickening
        # There are two solutions, one thin and one thick, so we get the thin case w/ the thicken condition
        if (abs(stress_tot - yield_str)/yield_str <= 0.001) and thicken:
            break

    return wall_t
