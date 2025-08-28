import pandas as pd
import numpy as np
import scipy.optimize as sci
import math

class JetA:
    M = .166 # kg/mol
    R = 8.3145 / M # J/kgK

    # https://pubs.acs.org/doi/pdf/10.1021/ie00040a045
    # Used their analytically estimated Jet-A critical pressure
    # Averaged their experimentally determined Jet-A critical temperature
    p_crit = 2_340_608 # Pa
    temp_crit = 662.039 # K

    # https://ntrs.nasa.gov/api/citations/19890007646/downloads/19890007646.pdf
    aGas1 = np.array([0.20869217e1, .13314965, -.81157452e-4, .29409286e-7, -.65195213e-11]) # for temp <= 1000K
    aGas2 = np.array([.24880201e2, .78250048e-1, -.31550973e-4, .578789e-8, -.39827968e-12]) # for temp > 1000K
    aLiquid = np.array([.19049613e2, -.16918532e-1, .63022035e-3, -.13336577e-5, .94335638e-9])

    # https://github.com/gpavanb-old/GroupContribution/blob/master/data/Init_Data/posf10325_surr_init.xlsx
    mole_fractions_dict = {"124 tmb" : .15, "iso-dodecane" : .3, "n-undecane" : 0.2,"pentyl-cyclohexane" : .35}

    acentric = -np.inf # defined & called in JetA.get_acentric_factor() ; DO NOT directly access this field

    # https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir6659.pdf page 59
    @staticmethod
    def get_rho_l(temp):
        # data = pd.read_excel('Jet-A Compressed Liquid Densities.xlsx',usecols=[1,2,3],skiprows=1)

        rhoSlope = (0.685-0.835)/(470-270) # From old data. Would be nice to update
        return (rhoSlope * (temp - 270) + 0.835) * 1000 #kg/m^3
    
    @staticmethod
    def get_kinematic_viscosity(temp):
        # Kinematic viscosity
        viscK_c =  1.7825 * (1700 * math.exp(-0.026*(temp)) + 0.27) / 1000000 # Done with eyeballed curve fit. Should be updated to use a more rigorous method
        return viscK_c

    @staticmethod
    def get_conductivity(temp):
        # Conductivity
        # Temperature in Celcius
        # Cond_300K = 0.1150 W/mK, Cond_550K = 0.076 W/(m*K), Pg. 82
        condSlope = (0.076-0.1150)/(550-300)

        # return condSlope * (temp - 300) + 0.1150 #W/(m*K)
        return .6
      
    @staticmethod
    def get_rho_v(temp, pressure):
        return pressure / (JetA.R * temp * JetA.get_Z(temp,pressure)) # state equation
    
    # https://asmedigitalcollection.asme.org/energyresources/article/136/1/012903/366580/Development-of-a-New-Correlation-of-Gas, equation 39
    @staticmethod
    def get_Z(temp, pressure):
        T_r = temp/JetA.temp_crit
        p_r = pressure/JetA.p_crit
        term1 = (.702*np.exp(-2.5*T_r))*p_r**2
        term2 = (5.524*np.exp(-2.5*T_r))*p_r
        term3 = (.044*T_r**2-.164*T_r+1.15)
        Z = term1 - term2 + term3
        if Z < .8 or Z > 1.2: # sanity check
            Z = 1
        return Z

    # a measure of how nonspherical molecules are
    # https://aiche.onlinelibrary.wiley.com/doi/epdf/10.1002/aic.690300615
    # lazily initialized (only runs the first time, since acentric factor does not depend on pressure, temperature, etc)
    @staticmethod
    def get_acentric_factor():
        M = JetA.M*1000 # kg/mol to g/mol

        if JetA.acentric > 0:
            return JetA.acentric
        
        M_arr = np.array([1, M, M**2, M**3, 1/M]) # used in equation 1

        # coefficients
        C_tb = np.array([1.33832e2, 3.11349, -7.08978e-3, 7.69085e-6, -1.12731e3])
        C_S = np.array([6.6405e-1, 1.4813e-3, -5.07021e-6, 6.21414e-9, -8.45218])
        C_omega = np.array([-1.56752*10,1.22751,9.96848e-3, -2.04742e-5, -6.90883*10])
        a_omega = np.array([-1.16044e3, 3.48210, 2.78317e4, -2.05257e2, 4.55767e-1, -7.13722e4, 5.08888e2, -6.10273e-1, -1.68712e-3])
        b_omega = np.array([1.89761, 2.41662e-2, -2.67462e2, 2.06071, -5.22105e-3, 7.6607e2, -5.75141, 8.66667e-3, 1.75189e-5])
        
        # normal boiling temp deviations
        Tb = 462.039 # normal boiling temp (K)
        Tba = np.dot(C_tb, M_arr) # normal boiling temp, K
        deltaTb = Tba - Tb
        
        # specific gravity deviations
        S = JetA.get_rho_l(20+273.15)/1000 # liquid density at 20 degrees C (and  1 atm) over water density at 5 degrees C
        Sa = np.dot(C_S, M_arr)
        deltaS = S-Sa

        # solve for acentric factor
        theta_A = np.dot(C_omega, M_arr)
        A_arr = a_omega+b_omega*JetA.M # equation 5
        A_arr = np.insert(A_arr, 0, 1)
        deltaST_arr = np.array([1, deltaS, deltaTb, deltaS**2, deltaS*deltaTb, deltaTb**2, deltaS**3, deltaS**2*deltaTb, deltaS*deltaTb**2,deltaTb**3])
        JetA.acentric = (theta_A + np.dot(A_arr, deltaST_arr))/JetA.temp_crit

        return JetA.acentric
    
    @staticmethod
    def get_surface_tension():
        return .1 # water value

    @staticmethod
    def get_cp_vapor(temp):
        tempVec = np.array([0, temp, temp**2, temp**3, temp**4])
        if temp <= 1000:
            return JetA.R*np.dot(JetA.aGas1, tempVec)
        else:
            return JetA.R*np.dot(JetA.aGas2, tempVec)

    @staticmethod
    def get_c_liquid(temp):
        tempVec = np.array([0, temp, temp**2, temp**3, temp**4])
        return JetA.R*np.dot(JetA.aLiquid, tempVec)
    
    # fix for temp dependency later
    # if can find hfg at a specific temp, then can integrate cpdT for liquid and gas phases from that temp to entered temp
    # then subtract those & add the known hfg value to that
    # https://web.stanford.edu/group/haiwanglab/HyChem/approach/Report_Jet_Fuel_Thermochemical_Properties_v6.pdf
    @staticmethod
    def get_h_fg():
        return .36e6
    
    # https://www.sciencedirect.com/science/article/pii/S0017931016309887?via%3Dihub#s0100, B4
    @staticmethod
    def get_v_dynamic_viscosity(temp):
        # zeta = .176*(JetA.temp_crit/(JetA.M**3*JetA.p_crit**4))**(1/6)*10**(-10/3)
        # temp_r = temp / JetA.temp_crit # reduced temperature
        # print(f"reduced temp: {temp_r}")
        # return (.807*temp_r**.618-.357*np.exp(-.449*temp_r)+.34*np.exp(-4.058*temp_r)+.018) / zeta
        return 1.5e-5 # dummy value ; this is water viscosity at 2,000K (actual equations aren't working)
    
    @staticmethod
    def get_l_dynamic_viscosity(temp):
        # tmb: https://dhc-solvent.de/downloads/DHC_SDS_042_en.pdf
        tmb_slope = (.528 - .727) / 30_000
        tmb_intercept = .727 * 1e-3 - 30 * tmb_slope

        # iso-dodecane: from unreliable sources ; an approximate value

        # Pentyl-cyclohexane : https://materials.springer.com/interactive?systemId=20566&propertyId=Dynamic%20Viscosity
        pentyl_dict = {383.15: 510, 378.15:540,373.15:570,368.15:600,363.15:630,
                       358.15:670,353.15:710,348.15:750,343.15:790,338.15:843,333.15:898,
                       328.15:958,323.15:1026,318.15:1103,313.15:1191,
                       308.15:1295,303.15:1418,298.15:1560,293.15:1723,288.15:1910,
                       283.15:2120,278.15:2360,273.2:2640,268.15:2690,263.15:3330}
        pentyl = JetA.interpolate(pentyl_dict, temp)/1000

        # n-undecane : https://pubchem.ncbi.nlm.nih.gov/compound/Undecane#section=Stability-Shelf-Life
        viscosities = {"124 tmb": tmb_intercept + temp*tmb_slope, "iso-dodecane": .001, "n-undecane":.001098,
                       "pentyl-cyclohexane" : pentyl}
        viscosity = 0
        for key in JetA.mole_fractions_dict.keys():
            viscosity += viscosities[key] * JetA.mole_fractions_dict[key]
        # to do later : add 'lucas' correction
        return viscosity

    # https://aiche.onlinelibrary.wiley.com/doi/epdf/10.1002/aic.690210313
    @staticmethod
    def get_saturation_temp(pressure):
        PR = pressure / JetA.p_crit # reduced pressure
        omega = JetA.get_acentric_factor()
        saturation_temp = sci.fsolve(func=JetA._saturation_temp_eqn, x0=[.1], args=(PR, omega))[0] * JetA.temp_crit
        return saturation_temp
        
    # numerical solver should make this equal 0
    @staticmethod
    def _saturation_temp_eqn(TR, PR, omega):
        f0 = 5.92714 - 6.09648 / TR - 1.28862 * np.log(TR) + .169347 * TR**6
        f1 = 15.2518 - 15.6875 / TR - 13.4721 * np.log(TR) + .43577 * TR**6
        return (f0 + omega * f1) - np.log(PR)

    # interpolates the table, represented as a dictionary of values
    @staticmethod
    def interpolate(dict, value):
        maxBelow = -np.inf
        minAbove = np.inf
        for key in dict.keys():
            if key == value:
                return dict[key]
            elif key > value:
                if key < minAbove:
                    minAbove = key
            else:
                if key > maxBelow:
                    maxBelow = key
        if maxBelow != -np.inf and minAbove != np.inf:
            slope = (dict[minAbove] - dict[maxBelow]) / (minAbove-maxBelow)
            return dict[maxBelow]+slope*(value-maxBelow)
        elif minAbove == np.inf:
            largestKey = max(dict.keys())
            largestKeyValue = dict[largestKey]
            newDict = dict.copy()
            del newDict[largestKey]
            nextLargestKey = max(newDict.keys())
            slope = (largestKeyValue-newDict[nextLargestKey])/(largestKey-nextLargestKey)
            difference = value-largestKey
            print(f"Warning : Extrapolateding up by {difference}")
            return largestKeyValue + slope * difference
        else:
            smallestKey = min(dict.keys())
            smallestKeyValue = dict[smallestKey]
            newDict = dict.copy()
            del newDict[smallestKey]
            nextSmallestKey = min(newDict.keys())
            slope = (newDict[nextSmallestKey]-smallestKeyValue)/(nextSmallestKey-smallestKey)
            difference = smallestKey - value
            print(f"Warning: extrapolating down by {difference}")
            return smallestKeyValue - slope * difference