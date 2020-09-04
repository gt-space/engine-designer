import math
from .findPrinciples import findPrinciples

def fos(R, t, meanT, dT):
    # Constants
    # si units
    Di = R / 2  # Diameter of engine station
    Do = Di + 2 * t # Outer diameter
    ri = Di / 2 # Inner radius
    ro = Do / 2 # Outer radius
    r_max = ri # Max radius
    Patm = 101325 # Atmospheric pressure in Pa
    Pi = 1.793e+6 # Gas pressure (Pa)
    Po = 2.344e+6 # Coolant pressure (Pa) P_waterflow = 137895 Pa; P_fuel_run = 2.344e+6
    F_nozzle = -3336.17 # Thrust

    # Thermal Stresses
    # Props from http://www-ferp.ucsd.edu/LIB/PROPS/PANOS/cu.html
    alpha = (13.251 + 6.903e-3 * meanT + 8.53063e-7 * meanT) * 10**-6 # CTE Currenty taking mean value
    E = 129.8e9 # al 6061 E = 68.9e9 Pa, al 7075 E = 71.7e9 Pa
    v = 0.34 # Possions ratio
    sigma_thermal = (2 * alpha * E * dT/(1 - v)) / 1000000

    kd_factor = 103 / 276 # Tmax = 204C
    Sy = 1570.3 - 14.184 * meanT + 5.641e-2 * meanT**2 - 1.0592e-4 * meanT**3 + 9.2881e-8 * meanT**4 - 3.086e-11 * meanT**5 # MPA Cu Annealed yield strength (MPa)
    Su = 191.31 + 0.65634 * meanT -1.85e-3 * meanT**2 + 1.0185e-6 * meanT**3 # MPA Cu Annealed ultimate strength (MPa)
    fos_desired = 1.1 # Factor of safety desired
    Kt = 1
    k_cond = 330 # k_cond_Cu110 = 330

    # STRESS CALCULATIONS

    # thick-walled vessel
    def sigma_t(r):
        # Calculate the tension stress of the chamber wall (Pa)
         return (Pi * ri**2 - Po * ro**2 - (ri**2) * (ro**2) * (Po - Pi) / r**2) / (ro**2 - ri**2)
    def sigma_r(r):
        # Calculate the radial stress of the chamber wall (Pa)
         return (Pi * ri**2 - Po * ro**2 + (ri**2) * (ro**2) * (Po - Pi) / r**2) / (ro**2 - ri**2)
    sigma_l = F_nozzle / (math.pi * (ro**2 - ri**2))

    # thin-walled approximation
    def sigma_t_approx(t):
        # Hoop stress for a thin wall (D/t > 30) can be approximated with the following
        return (Pi - Po) * (Di + t)/(2 * t)

    sigma_hoop = sigma_t_approx(t) / 1000000 # Pa to MPa

    # Constants
    sigma_x = sigma_l / 1000000 # Pa to MPa
    sigma_y = sigma_t(r_max) / 1000000 + sigma_thermal # Pa to MPa
    sigma_z = sigma_r(r_max) / 1000000 # Pa to MPa
    tau_xy = 0
    tau_yz = 0
    tau_zx = 0

    # Cubic Function
    def f(sigma):
        return sigma**3 - (sigma_x + sigma_y + sigma_z) * sigma**2 + (sigma_x * sigma_y + sigma_x * sigma_z + sigma_y * sigma_z - tau_xy**2 - tau_yz**2 - tau_zx**2) * sigma - (sigma_x * sigma_y * sigma_z + 2 * tau_xy * tau_yz * tau_zx - sigma_x * tau_yz**2 - sigma_y * tau_zx**2 - sigma_z * tau_xy**2)

    # Find principle stresses (MPa)
    [sigma_1, sigma_2, sigma_3] = findPrinciples(sigma_x, sigma_y, sigma_z, tau_xy, tau_yz, tau_zx)

    # Principle shears ()
    tau_12 = (sigma_1 - sigma_2) / 2
    tau_23 = (sigma_2 - sigma_3) / 2
    tau_13 = (sigma_1 - sigma_3) / 8

    # Maximum Shear Stress Yield Criterion
    def ny_max_shear(Sy, s1, s3):
        return Sy/(s1 - s3)

    # Octahedral Shear Stress Yield Criterion (Distortion Energy)
    def ny_dist_eng(Sy, s1, s2, s3):
        return Sy / (((s1 - s2)**2 + (s2 - s3)**2 + (s3 - s1)**2) / 2)**0.5

    # Calculate factors of safety using both theories of failure
    ny_conserv = ny_max_shear(Sy, sigma_1, sigma_3)
    ny_precise = ny_dist_eng(Sy, sigma_1, sigma_2, sigma_3)

    nu_conserv = ny_max_shear(Su, sigma_1, sigma_3)
    nu_precise = ny_dist_eng(Su, sigma_1, sigma_2, sigma_3)

    return (ny_conserv, ny_precise, nu_conserv, nu_precise)
