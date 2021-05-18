# alt_design.py : Defines engine class and applies relevant equations to build
# and simulate and engine geometry given input parameters

# Uses a fixed exit area and then back solves for thrust
# NOTE: This code has not been updated to work with newest version of main.py
# Currently tabled until needed

# OUTPUTS:
# engineContour is a numpy array with the following column format:
# [R, Z]
# engineProps is a numpy array with the following column format:
# [R, Z, pip, aeat, mach, cf, ivac, isp, p, t, rho, h, u, mw, cp, gam, son,
# vis, cond, pran, cpfz, condfz, pranfz]

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj
from contour import get_contour
from CEA_properties import get_props
from constants import Constants

np.set_printoptions(threshold=sys.maxsize) # Print setting for debugging arrays

class AltEngine:
    # Upon declaration of a new engine:
    def __init__(self, exit_diam, P_inj, conRat, MR=2, divAng=15, conAng=35, radRat=0.7, LStar=1.05):
        self.exit_diam = exit_diam # Exit diameter [m]
        self.P_inj = P_inj # P_inj: Injector face pressure [Bar]
        self.conRat = conRat # conRat: Enigne contraciton ratio
        self.MR = MR # Mixture ratio by weight (ox/fuel)
        self.divAng = divAng # Divergence half angle [Deg]
        self.conAng = conAng # Convergence half angle [Deg]
        self.radRat = radRat # radRat: R/Rmax (BETWEEN 0 & 1)
        self.LStar = LStar # Characteristic length [m]

    # Main design function
    def design_engine(self):
        # Call in global constants
        constants = Constants()

        # ==== Specify Parameters ====
        P_inj_psi = self.P_inj * constants.psi_to_bar # Psi is used in many CEA funcitons
        P_rat = self.P_inj/constants.P_amb # Inverse pressure ratio of overall expansion
        # P_rat *= 2 # Correction factor for altitude (~5km ideal operating alt for 2)

        # Define RocketCEA object
        ispObj = CEA_Obj( oxName='LOX', fuelName='JetA', fac_CR=self.conRat)

        # ==== Get Expansion Ratio ====

        # Throat pressure ratio (injector pressure / throat pressure)
        pip_t = ispObj.get_Throat_PcOvPe(P_inj_psi, self.MR)
        P_t = (self.P_inj / pip_t) * constants.bar

        # Get the area ratio where optimal expansion happens
        exp_rat = ispObj.get_eps_at_PcOvPe(P_inj_psi, self.MR, P_rat)

        # ==== Get Throat Radius ====

        A_t = math.pi * (self.exit_diam ** 2) / 4 / exp_rat

        R_t = math.sqrt(A_t / math.pi) # Calculate Throat Radius (m)

        # Get Mass Flow Rate

        T = ispObj.get_Temperatures(P_inj_psi, self.MR, exp_rat, frozen=0, frozenAtThroat=0) # Get temperature array

        T_t = (T[1]-32)*5/9 + constants.C_to_K # Throat temperature (F to K)
        T_c = (T[0]-32)*5/9 + constants.C_to_K # Throat temperature (F to K)

        (MW, gam) = ispObj.get_Throat_MolWt_gamma(P_inj_psi, self.MR) # Get molecular weight and specific heat ratio
        R_spec = constants.Ru / MW

        self.mDot_tot = (A_t * self.P_inj * constants.bar / math.sqrt(T_c)) * math.sqrt(gam/R_spec) * ((gam + 1) / 2)**(-(gam + 1) / (2 * (gam - 1))) # Mass flow rate assuming choked flow at throat (kg/s)
        self.mDot_o = (self.mDot_tot/(self.MR + 1)) * self.MR # Oxidizer mass flow (kg/s)
        self.mDot_f = self.mDot_tot - self.mDot_o # Fuel mass flow (kg/s)

        # Check Calcs:
            # k = gam
            #
            # self.mDot_tot_hh = A_t * P_t * k* math.sqrt((2 / (k+1))**((k+1)/(k-1)))/math.sqrt(k * R_spec * T_t)
            # print(self.mDot_tot_hh)
            #
            # self.mDot_tot_other = A_t * P_t * math.sqrt(gam/(R_spec * T_t)) # Old mDot equation, should be similar
            #
            # print(self.mDot_tot_other)

        M = ispObj.get_MachNumber(P_inj_psi, self.MR, exp_rat) # Calculate mach number at exit

        v_sonic = ispObj.get_SonicVelocities(P_inj_psi, self.MR, exp_rat) # Calculate sonic velocity [ft/sec]

        v_sonic = v_sonic[2]*constants.ft_to_m # [ft/sec to m/s]


        # Check expansion ratio:
        k = gam

        E = (1/M) * ((1+(k-1)*(M**2)/2)/((k+1)/2)) ** ((k+1)/(2*(k-1)))
        print(E, exp_rat)

        correction_factor = (1 + math.cos(self.divAng*math.pi/180))/2 # Correction factor for nozzle exit velocity, using conic appx
        v_exit = v_sonic * M # Sonic Vel * Mach @ Exit [m/s]

        self.thrust = self.mDot_tot * (correction_factor * v_exit)  # Calculate total mass flow based on V_exit [kg/s]


        self.conAng *= math.pi / 180 # Nozzle contraction angle [degrees to rad]

        # Generate engine contour from external function:
        (self.engineContour, self.chBarrel, self.nozzleContour, R_tCurve, self.throatInd, self.conLeadInRadius) = getContour(R_t, self.LStar, self.conRat, self.conAng, self.divAng, self.radRat, exp_rat)

        # Generate engine property array from helper function
        # This creates a 200x23 array describing various properties along the nozzle.
        # See getProperties.py for details
        self.engineProps = getProps(self.chBarrel, self.nozzleContour, self.throatInd, ispObj, P_inj_psi, self.MR, A_t)

        # Check for bad design (This occurs when thrust is way higher than chamber pressure should be and barrel becomes negative)
        if self.chBarrel[1,1] < 0:
            raise Exception("Thrust is too high for given chamber pressure to create a geometry with this method. \
                Try lowering thrust, increasing chamber pressure and/or decreasing contraction ratio.")

        # Performance parameters (not used in calculations)
        self.C_f = ispObj.get_PambCf(constants.psi_to_atm, P_inj_psi, self.MR, exp_rat)
        self.C_star = ispObj.get_Cstar(P_inj_psi, self.MR) * constants.ft_to_m * 0.975 # Calculate C*, m/s, apply 0.975 correction factor per H&H (pp. 70)
