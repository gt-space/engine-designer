# ----==== GENERAL INFO ====----
# Last Updated: 7/24/2020
# This script was adapted from an earlier MATLAB version. It is designed
# to size an engine based on a few parameters and return a data array
# describing various conditions at points along this created engine. It does
# this by generating a contour vector for the engine and then querying results
# from the RocketCEA library. To run this, you will need the following:
#   * Python 3 installed on your machine (Script was developed and tested with 3.7.3)
#   * RocketCEA installed
#   * numpy installed
#   * gofortran, which is a dependancy for rocketcea
#   * matplotlib if you want to see or make plots
# Installation Guide: https://rocketcea.readthedocs.io/en/latest/quickstart.html
#
# NOTE: RocketCEA website documentation references combustion end pressure in its
# funcitons, but it should say injector face pressure.

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj
from engine_designer.getContour import getContour
from engine_designer.getProperties import getProps

np.set_printoptions(threshold=sys.maxsize) #Print setting for debugging

# MR: ox/fuel
# thrust_nom: Design thrust (N)
# P_inj: Injector face pressure (Bar)
# conRat: Enigne contraciton ratio
# divAng: Divergence half angle (Deg)
def design_engine(thrust_nom, P_inj, conRat, MR = 1.8, divAng = 15):
    # ==== Propellants: Kerosene and Liquid Oxygen ====
    LStar = 1.05 # Characteristic length (m)
    # ==== Define Constants ====
    g0 = 9.81 # Gravity (m/s^2)
    bar = 100000 # 1 Bar in Pa
    Ru = 8314.46 # Universal Gas Constant (J/kmolK)

    # ==== Specify Parameters ====
    P_inj_psi = P_inj * 14.5038 # Psi is used in many CEA funcitons
    P_amb = 1.01325 # Ambient Pressure (Bar)
    P_rat1 = P_amb/P_inj # Pressure ratio of overall expansion
    P_rat2 = P_inj/P_amb # Inverse pressure ratio of overall expansion


    ispObj = CEA_Obj( oxName='LOX', fuelName='JetA', fac_CR=conRat)

    pip_t = ispObj.get_Throat_PcOvPe(P_inj_psi, MR)

    # Find some Parameters...
    # get the area ratio where optimal expansion happens
    expRat = ispObj.get_eps_at_PcOvPe(P_inj_psi, MR, P_rat2)

    C_f = ispObj.get_PambCf(14.69594878, P_inj_psi, MR, expRat)

    M = ispObj.get_MachNumber(P_inj_psi, MR, expRat)

    sonic_v = ispObj.get_SonicVelocities(P_inj_psi, MR, expRat)
    sonic_v = sonic_v[2]*0.3048

    # ====Size & Calculate Engine Parameters====
    # This creates a 200x22 array describing various properties along the nozzle.
    # See  getProperties.py for details

    nozCorrFactor = (1 + math.cos(divAng*math.pi/180))/2 # Correction factor for nozzle exit velocity, using conic appx
    V_exit = sonic_v * M # Sonic Vel * Mach @ Exit (m/s)
    Cf_exit = C_f # Exit thrust coeficient
    mDot_tot = thrust_nom/(nozCorrFactor * V_exit) # Calculate total mass flow based on V_exit (kg/s)
    mDot_o = (mDot_tot/(MR + 1)) * MR # Oxidizer mass flow (kg/s)
    mDot_f = mDot_tot - mDot_o # Fuel mass flow (kg/s)

    T = ispObj.get_Temperatures(P_inj_psi, MR, expRat, frozen=0, frozenAtThroat=0)

    T_c = (T[0]-32)*5/9 # Chamber temperature (F to C)
    T_t = (T[1]-32)*5/9 # Throat temperature (F to C)

    (MW, gam) = ispObj.get_Throat_MolWt_gamma(P_inj_psi, MR) # Get molecular weight and specific heat ratio

    A_t = (mDot_tot/((P_inj / pip_t) * bar)) * math.sqrt((Ru * T_t)/(MW * gam)) # Area at throat (m^2)
    C_star = ispObj.get_Cstar(P_inj_psi, MR) * 0.3048 * 0.975 # Calculate C*, m/s, apply 0.975 correction factor per H&H (pp. 70)
    R_t = math.sqrt(A_t/math.pi) # Calculate Throat Radius (m)
    R_c = R_t*math.sqrt(conRat) # Chamber Radius (m)
    radRat = 0.7 # R2/R2max - Refer to RPA User Manual : BETWEEN 0 & 1
    conAng = 35 * math.pi / 180 # Nozzle contraction angle (degrees to rad)

    # Generate engine contour from external function:
    (engineContour, chBarrel, nozzleContour, R_tCurve, throatInd) = getContour(R_t, LStar, conRat, conAng, divAng, radRat, expRat)

    # Generate engine property array from helper function
    engineProps = getProps(chBarrel, nozzleContour, throatInd, ispObj, P_inj_psi, MR, A_t)

    return (engineContour, engineProps)
