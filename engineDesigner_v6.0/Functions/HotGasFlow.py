from scipy.optimize import brentq
import numpy as np

def area_mach_relation(M, gamma):
    return (1/M)*((2/(gamma+1))*(1+(gamma-1)/2*M**2))**((gamma+1)/(2*(gamma-1)))

def solve_mach(area_ratio, gamma, branch):
    # Throat
    if np.isclose(area_ratio, 1.0, atol=1e-8):
        return 1.0
    
    def residual(M):
        return area_mach_relation(M, gamma) - area_ratio

    if branch == "subsonic":
        return brentq(residual, 1e-6, 0.999999)

    elif branch == "supersonic":
        return brentq(residual, 1.000001, 20.0)

#Using isentropic + Mach relations 
def static_temperature(T0, gamma, M):
    return T0 / (1 + (gamma - 1)/2 * M**2)


def adiabatic_wall_temperature(T, T0, Pr):
    recovery_factor = Pr**(1/3)
    return T + recovery_factor * (T0 - T)


def static_pressure(Pc, gamma, M):
    return Pc * (1 + (gamma - 1)/2 * M**2)**(-gamma/(gamma - 1))


def density(P, R, T):
    return P / (R * T)


def gas_velocity(M, gamma, R, T):
    return M * np.sqrt(gamma * R * T)

    