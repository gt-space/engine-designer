import numpy as np
from Geometry.Contour import getCurrentContour
from Geometry.Contour import getCurrentEngineParameters
from EngineParameters import(fuel, oxidizer, Pc, MR, Cstar_eff)
from Functions import HotGasFlow
from rocketcea.cea_obj import CEA_Obj

# --- 1. CALCULATE GAS PROPERTIES USING CEA ---
def getGasPropertiesCEA():

    geometry_params = getCurrentEngineParameters()
    eps_exp = geometry_params["eps"]   

    rocket = CEA_Obj(oxName=oxidizer, fuelName=fuel)
    Isp, Cstar, T0, MW, gamma = rocket.get_IvacCstrTc_exitMwGam(Pc, MR, eps_exp, 1, 0)
    cp, mu, k, Pr= rocket.get_Chamber_Transport(Pc, MR, eps_exp, 1)

    #Calculate Cstar actual to be used in Bartz Later
    Cstar_actual = Cstar*Cstar_eff

    # Convert RocketCEA default units to SI
    T0 *= 5.0 / 9.0      # deg R -> K
    Cstar_actual *= 0.3048      # ft/s -> m/s
    mu *= 1e-4 #millipoise -> Pa*s
    k  *= 0.4184 #mcal/(cm*s*K) -> W/(m*K)
    cp *= 4184.0 #cal/(gK) -> J/(kgK)

    # Specific gas constant
    R = 8314.4626 / MW

    return {
        "T0": T0,
        "MW": MW,
        "gamma": gamma,
        "cp_g": cp,
        "mu_g": mu,
        "k_g": k,
        "Pr": Pr,
        "R": R,
        "Cstar_actual": Cstar_actual
    }

# --- 3. CALCULATE FLOW PROPERTIES ALONG ENGINE ---

def getHotGasFlow(contour):

    # Convert chamber pressure to SI
    Pc_Pa = Pc * 6894.75729

    #Get engine geometry
    A = contour["A"]
    throat_index = contour["throat_index"]
    N = contour["N"]

    #Get hot gas properties
    gas = getGasPropertiesCEA()
    gamma = gas["gamma"]
    T0 = gas["T0"]
    Pr = gas["Pr"]
    R = gas["R"]

    #Initialize Arrays With Outputs that will be calculated later
    area_ratio = np.zeros(N)
    M = np.zeros(N)
    T = np.zeros(N)
    Taw = np.zeros(N)
    P = np.zeros(N)
    rho = np.zeros(N)
    V = np.zeros(N)

    #1D Quasi Flow Analysis
    for i in range(N):
        
        #Solve for area ratio
        area_ratio[i] = A[i] / A[throat_index]
        
        #Solve for mach number using Area-Mach Relation
        if i < throat_index:
            M[i] = HotGasFlow.solve_mach(area_ratio[i],gamma,"subsonic")
        elif i == throat_index:
            M[i] = 1.0
        else:
            M[i] = HotGasFlow.solve_mach(area_ratio[i], gamma, "supersonic")
        
        #Solve for static temperature
        T[i] = HotGasFlow.static_temperature(T0, gamma, M[i])

        #Solve for adabiatic wall temperature
        Taw[i] = HotGasFlow.adiabatic_wall_temperature(T[i], T0, Pr)
        
        #Solve for static pressure
        P[i] = HotGasFlow.static_pressure(Pc_Pa, gamma, M[i])
        
        #Solve for density
        rho[i] = HotGasFlow.density(P[i], R, T[i])
        
        #Solve for Gas Velocity
        V[i] = HotGasFlow.gas_velocity(M[i], gamma, R, T[i])

    return {
        "M": M,
        "T": T,
        "Taw": Taw,
        "P": P,
        "rho": rho,
        "V": V
        }
