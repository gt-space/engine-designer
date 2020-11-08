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
from .getContour import getContour
from .getProperties import getProps

np.set_printoptions(threshold=sys.maxsize) #Print setting for debugging

# MR: ox/fuel
# thrust_nom: Design thrust (N)
# P_inj: Injector face pressure (Bar)
# conRat: Enigne contraciton ratio
# divAng: Divergence half angle (Deg)
# Characteristic length (m)
# radRat: R2/R2max - Refer to RPA User Manual : BETWEEN 0 & 1
class Engine:
    def __init__(self, thrust_nom, P_inj, conRat, MR=2, divAng=15, conAng=35, radRat=0.7, LStar=1.05):
        self.thrust_nom = thrust_nom
        self.P_inj = P_inj
        self.conRat = conRat
        self.MR = MR
        self.divAng = divAng
        self.conAng = conAng
        self.radRat = radRat
        self.LStar = LStar

    def design_engine(self):
        # ==== Propellants: Kerosene and Liquid Oxygen ====
        # ==== Define Constants ====
        self.g0 = 9.81 # Gravity (m/s^2)
        self.bar = 100000 # 1 Bar in Pa
        self.Ru = 8314.46 # Universal Gas Constant (J/kmolK)

        # ==== Specify Parameters ====
        self.P_inj_psi = self.P_inj * 14.5038 # Psi is used in many CEA funcitons
        self.P_amb = 1.01325 # Ambient Pressure (Bar)
        self.P_rat = self.P_inj/self.P_amb # Inverse pressure ratio of overall expansion

        # Define RocketCEA object
        ispObj = CEA_Obj( oxName='LOX', fuelName='JetA', fac_CR=self.conRat)

        # ==== Get Expansion Ratio ====

        # Throat pressure ratio (injector pressure / throat pressure)
        self.pip_t = ispObj.get_Throat_PcOvPe(self.P_inj_psi, self.MR)

        # Get the area ratio where optimal expansion happens
        self.expRat = ispObj.get_eps_at_PcOvPe(self.P_inj_psi, self.MR, self.P_rat)

        # ==== Get Throat Radius ====

        self.M = ispObj.get_MachNumber(self.P_inj_psi, self.MR, self.expRat)

        sonic_v = ispObj.get_SonicVelocities(self.P_inj_psi, self.MR, self.expRat)
        self.sonic_v = sonic_v[2]*0.3048

        self.nozCorrFactor = (1 + math.cos(self.divAng*math.pi/180))/2 # Correction factor for nozzle exit velocity, using conic appx
        self.V_exit = self.sonic_v * self.M # Sonic Vel * Mach @ Exit (m/s)

        self.mDot_tot = self.thrust_nom/(self.nozCorrFactor * self.V_exit) # Calculate total mass flow based on V_exit (kg/s)
        self.mDot_o = (self.mDot_tot/(self.MR + 1)) * self.MR # Oxidizer mass flow (kg/s)
        self.mDot_f = self.mDot_tot - self.mDot_o # Fuel mass flow (kg/s)

        self.T = ispObj.get_Temperatures(self.P_inj_psi, self.MR, self.expRat, frozen=0, frozenAtThroat=0) # Get temperature array

        self.T_t = (self.T[1]-32)*5/9 + 273.15 # Throat temperature (F to K)

        (self.MW, self.gam) = ispObj.get_Throat_MolWt_gamma(self.P_inj_psi, self.MR) # Get molecular weight and specific heat ratio

        self.A_t = (self.mDot_tot/((self.P_inj / self.pip_t) * self.bar)) * math.sqrt((self.Ru * self.T_t)/(self.MW * self.gam)) # Area at throat (m^2)

        self.R_t = math.sqrt(self.A_t / math.pi) # Calculate Throat Radius (m)

        self.conAng = self.conAng * math.pi / 180 # Nozzle contraction angle (degrees to rad)

        # Generate engine contour from external function:
        (self.engineContour, self.chBarrel, self.nozzleContour, self.R_tCurve, self.throatInd) = getContour(self.R_t, self.LStar, self.conRat, self.conAng, self.divAng, self.radRat, self.expRat)

        # Generate engine property array from helper function
        # This creates a 200x22 array describing various properties along the nozzle.
        # See getProperties.py for details
        self.engineProps = getProps(self.chBarrel, self.nozzleContour, self.throatInd, ispObj, self.P_inj_psi, self.MR, self.A_t)

        # Check for bad design (This occurs when thrust is way higher than chamber pressure should be and barrel becomes negative)
        if self.chBarrel[1,1] < 0:
            raise Exception("Thrust is too high for given chamber pressure to create a geometry with this method. \
                Try lowering thrust, increasing chamber pressure and/or decreasing contraction ratio.")

        # Performance parameters (not used in calculations)
        self.C_f = ispObj.get_PambCf(14.69594878, self.P_inj_psi, self.MR, self.expRat)
        self.C_star = ispObj.get_Cstar(self.P_inj_psi, self.MR) * 0.3048 * 0.975 # Calculate C*, m/s, apply 0.975 correction factor per H&H (pp. 70)
