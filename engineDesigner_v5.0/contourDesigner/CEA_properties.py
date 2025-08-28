# ==== Calculate Properties ====
# The array generated shall have the following header format:
# R, Z, pip, aeat, mach, cf, ivac, isp, p, t, rho, h, u, mw, cp, gam, son,
#   vis, cond, pran, cpfz, condfz, pranfz
# A description of each of these parameters can be found at:
# https://cearun.grc.nasa.gov/cgi-bin/CEARUN/plotParameters.cgi

import math
import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj


# Helper function for parsing the large output string. It might be possible to
# optimize all the loops, but for now it seems to work well enough.
def getValue(output, str, i=1, n=1):
    # i is the col you're trying to get from the row of string (first item is 1)
    # n the "nth" appearance of this value in the output string that you want.
    ind = output.find(str) + len(str)
    # l is the length of the value we're finding
    l = 1
    if str == "Cp, KJ/(KG)(K)":
        n += 1
    # rho is formatted oddly in the output string, so it needs to be handled
    # differently. This can probably be made neater, but for now it works.
    if str == "RHO, KG/CU M":
        while i > 1:
            # When we're not yet at the point
            while output[ind] == " ":
                ind += 1
            # passing some earlier numbers
            while output[ind] != " ":
                ind += 1
            if output[ind:ind+2] == " 0":
                ind += 2
            i -= 1
        # figure out and return value if you're at its index:
        exp = 0
        while output[ind] == " ":
            ind += 1
        while output[ind + l] != " ":
            l += 1
            if output[ind + l] == "+":
                exp = int(output[ind + l + 1])
                break
            elif output[ind + l] ==  "-":
                exp = -1 * int(output[ind + l + 1])
                break
        return float(output[ind:ind+l]) * (10 ** exp)
    while ind >= 0 and n > 1:
        ind = output.find(str, ind + len(str)) + len(str)
        n -= 1
    # Some numbers have exponents, so that must be taken into account
    while i > 1:
        # When we're not yet at the point
        while output[ind] == " ":
            ind += 1
        # passing some earlier numbers
        while output[ind] != " ":
            ind += 1
        i -= 1
    # figure out and return value if you're at its index:
    while output[ind] == " ":
        ind += 1
    while output[ind + l] != " ":
        l += 1
    return float(output[ind:ind+l])

def get_props(chBarrel, nozzleContour, throatInd, ispObj, P_inj_psi, MR, A_t, cstar_eff = 1):
    # Find the Chamber Barrel Properties. Idealized as linear from inj to comb end
    # Chamber barrel length
    numPTS = len(chBarrel)
    # full CEA output (used when other CEA methods are not available to gather specific data)
    output = ispObj.get_full_cea_output(P_inj_psi, MR, short_output=1, output="siunits")
    print(output)
    error

    #print(output)
    # Define chamber props
    chamberProps = chBarrel

    # These are strings to index the large output string. Any spaces you see are needed.
    dataItems = ["Pinj/P ", "Ae/At", "MACH NUMBER", "CF", "Ivac, M/SEC", "Isp, M/SEC", "P, BAR", "T, K",
    "RHO, KG/CU M", "H, KJ/KG", "U, KJ/KG", "M, (1/n)", "Cp, KJ/(KG)(K)", "GAMMAs", "SON VEL,M/SEC",
    "VISC,MILLIPOISE", "CONDUCTIVITY  ", "PRANDTL NUMBER"]

    # Items that do not have injector face data in the output string
    no_inj = ["Ae/At", "CF", "Ivac, M/SEC", "Isp, M/SEC"]

    # Items that appear more than once
    repeatedItems = ["Cp, KJ/(KG)(K)", "CONDUCTIVITY  ", "PRANDTL NUMBER"]

    # Generate barrel properties
    for i in range(len(dataItems)):
        if dataItems[i] == "Ae/At":
            inj = getValue(output, dataItems[i])
            end = getValue(output, dataItems[i])
        elif dataItems[i] in no_inj:
            inj =  0
            end = getValue(output, dataItems[i])
        else:
            inj = getValue(output, dataItems[i])
            end = getValue(output, dataItems[i], 2)
        chamberProps = np.column_stack((chamberProps, np.linspace(inj, end, numPTS)))
    for i in range(len(repeatedItems)):
        # We want to ensure we're getting the frozen values
        inj = getValue(output, repeatedItems[i], 1, 2)
        end = getValue(output, repeatedItems[i], 2, 2)
        chamberProps = np.column_stack((chamberProps, np.linspace(inj, end, numPTS)))


    nozzleProps = np.column_stack((nozzleContour, np.zeros((len(nozzleContour), len(dataItems) + len(repeatedItems)))))

    #generate converging nozzle properties. This can probably be optimized later.
    for a in range(throatInd):
        aRat = ((nozzleContour[a, 0] ** 2) * math.pi)/ A_t
        output = ispObj.get_full_cea_output(P_inj_psi, MR, subar=aRat, short_output=1, output="siunits")
        # Grab properties
        for i in range(len(dataItems)):
            if dataItems[i] in no_inj:
                prop = getValue(output, dataItems[i], 3)
            else:
                prop = getValue(output, dataItems[i], 4)
            nozzleProps[a, i+2] = prop
        for i in range(len(repeatedItems)):
            prop = getValue(output, repeatedItems[i], 4, 2)
            nozzleProps[a, i + 2 + len(dataItems)] = prop

    # generate diverging nozzle properties
    for a in range(throatInd,len(nozzleProps)):
        aRat = ((nozzleContour[a, 0] ** 2) * math.pi )/ A_t
        output = ispObj.get_full_cea_output(P_inj_psi, MR, eps=aRat, short_output=1, output="siunits")
        # Grab properties
        for i in range(len(dataItems)):
            if dataItems[i] in no_inj:
                prop = getValue(output, dataItems[i], 3)
            else:
                prop = getValue(output, dataItems[i], 4)
            nozzleProps[a, i+2] = prop
        for i in range(len(repeatedItems)):
            prop = getValue(output, repeatedItems[i], 4, 2)
            nozzleProps[a, i + 2 + len(dataItems)] = prop

    engineProps = np.concatenate((chamberProps, nozzleProps))

    return engineProps


# functions to convert from standard metric (kg, seconds, meters, Kelvin, Pascals, Joules) to units used in CEA or vice versa

# cea units: ft/s for cstar and velocities, psia for pressure, degR for temp, BTU/lbm for enthalpy,
# lbm/ft**3 for density, btu/lbm degR for specific energy capacity, millipoise for viscosity, mcal/cm-K-s for thermal conductivity

# metric units used : m/s for cstar and velocities, pascals for pressure, K for temp, J/kg for enthalpy,
# kg/m**3 for density, J/kgK for specific energy capacity, Pa*s for viscosity, J/m*K*s for thermal conductivity

def ceaToSI(number, value):
    cea_to_si_conversions = {"velocity": .3048, "pressure": 6894.76, "temperature": 5/9, 
                            "enthalpy": 2326, "density": 16.018, "specific heat": 2326 / (5 / 9),
                            "viscosity": 1e-4, "conductivity": 418.4}
    try:
        return number * cea_to_si_conversions[value]
    except KeyError:
        print("Unrecognized value string: " + value)
        return number

def siToCEA(number, value):
    cea_to_si_conversions = {"velocity": .3048, "pressure": 6894.76, "temperature": 5/9, 
                            "enthalpy": 2326, "density": 16.018, "specific heat": 2326 / (5 / 9),
                            "viscosity": 1e-4, "conductivity": .4184}
    try:
        return number / cea_to_si_conversions[value]
    except KeyError:
        print("Unrecognized value string: " + value)
        return number
