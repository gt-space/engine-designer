import math

# We want to find the maximum strain on the part. Assumption that this occurs at
# Midpoint of gas side wall.

# Strain will be the sum of these three components (in decreasing significance)
# 1. Thermal strain due to restricted expansion of the piece as a whole using average wall temp
# 2. Thermal strain due to relative temperature of the two sides of the wall
# 3. Strain from a pressure difference

def get_strain(R, t, meanT, dT, channel_w):
    strain_max = 0.3
    # Constants (si units)
    Patm = 101325 # Atmospheric pressure in Pa
    Pi = 1.793e+6 # Gas pressure (Pa)
    Po = 2.344e+6 # Coolant pressure (Pa) P_waterflow = 137895 Pa; P_fuel_run = 2.344e+6

    # Props from http://www-ferp.ucsd.edu/LIB/PROPS/PANOS/cu.html
    alpha = (13.251 + 6.903e-3 * meanT + 8.53063e-7 * meanT) * 10**-6 # CTE Currenty taking mean value
    E = 129.8e9 # al 6061 E = 68.9e9 Pa, al 7075 E = 71.7e9 Pa
    v = 0.34 # Possions ratio

    strain_t1 = -alpha * (meanT - 298)

    strain_t2 = -2 * alpha * dT/(1 - v)

    strain_p = 3*(Po - Pi)*(channel_w)**2/4*(t)**2 # Derived from max beam bending. Should double check calcs.

    strain = strain_t1 + strain_t2 + strain_p

    cycle_limit = -strain_max/strain # Extremely conservative cycle limit (likely capable of 4x cycles)
    # print("Strain: " + str(strain), "Cycle Limit:" + str(cycle_limit))

    return strain
