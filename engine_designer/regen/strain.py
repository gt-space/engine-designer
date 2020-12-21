import math

# We want to find the maximum strain on the part. Assumption is that this occurs at
# Midpoint of gas side wall.

# Strain will be the result of these three components
# 1. Thermal strain due to restricted expansion of the piece as a whole using average wall temp
# 2. Thermal strain due to relative temperature of the two sides of the wall
# 3. Strain from a pressure difference

# The thermal strain occurs due to a temperature increase which causes an intrinsic strain. This will
# bend the inner wall by theta_thermal. The pressure difference will induce a tensile stress in the material
# as it bows out into the chamber. The bend angle must be large enough to maintain static equilibrium, so if the
# derived minimum pressure bend angle is larger than the thermal bend angle, we will bend up to that pressure angle
# Otherwise the angle remains at the thermal one.

# Strain is calculated from this angle.

def get_strain(R, t, meanT, dT, channel_w, P_c, chamber_pressure):
    strain_max = 0.3
    # Constants (si units)
    Patm = 101325 # Atmospheric pressure in Pa
    Pi = chamber_pressure * 100000 # Gas pressure (Pa)
    Po = P_c * 100000 # Coolant pressure (Pa) P_waterflow = 137895 Pa; P_fuel_run = 2.344e+6

    # Props from http://www-ferp.ucsd.edu/LIB/PROPS/PANOS/cu.html
    alpha = (13.251 + 6.903e-3 * meanT + 8.53063e-7 * meanT) * 10**-6 # CTE Currenty taking mean value
    E = 129.8e9 # al 6061 E = 68.9e9 Pa, al 7075 E = 71.7e9 Pa
    E = 70e9 # Rough estimate for high temps
    v = 0.34 # Poisson's ratio for copper
    Sy = 1570.3 - 14.184 * meanT + 5.641e-2 * meanT**2 - 1.0592e-4 * meanT**3 + 9.2881e-8 * meanT**4 - 3.086e-11 * meanT**5 # MPA Cu Annealed yield strength (MPa)

    strain_t1 = alpha * (meanT - 298) # Wall strain due to temperature increase from ambient
    strain_t2 = 2 * alpha * dT/(1 - v) # Strain on inner side of wall due to temp difference across wall

    strain_thermal = strain_t1 + strain_t2 # Combined thermal strain

    # Bend angle caused by thermal strian:
    # theta_thermal = math.sqrt((-1/6 + math.sqrt((1/6)**2 - (7 * strain_thermal / 90)))/(7/180)) # Approximated theta with two term Taylor expansion: strain = 1/6 * theta^2 + 7/360 * theta^4

    # Old beam analysis method
    # I = t**3/12 # Adjusted bending moment (there is no b term as we are using pressure for the distributed load)
    # theta = (Po-Pi)*(channel_w)**2/(24*E*I)
    # strain_p = 3*(Po - Pi)*(channel_w)**2/4*(t)**2 # Derived from max beam bending

    # Retired Bending Calcs:
    # r_bend = t*Sy*1000000/(Po - Pi*t)
    # sigma_pressure = r_bend*(Po-Pi*t)/t - Pi
    # print(sigma_pressure/1000000, Sy)
    # theta_pressure = math.asin(channel_w/(2*r_bend))

    # Solve for minimum pressure bend angle
    theta_p = math.asin(channel_w*(Po - Pi)/(2*Sy*1000000*t))
    strain_p = theta_p/math.sin(theta_p) - 1

    # Determine total strain
    strain = strain_thermal + strain_p

    cycle_limit = strain_max/strain # Conservative cycle limit estimate (likely capable of 4x as many cycles)

    return strain, cycle_limit
