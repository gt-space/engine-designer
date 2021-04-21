# DESIGN.PY – Defines engine class and applies relevant equations to build
# and simulate and engine geometry given input parameters

# OUTPUTS
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
from contour import getContour
from CEA_properties import getProps
from constants import Constants

np.set_printoptions(threshold=sys.maxsize) # Print setting for debugging arrays

## DEFINE ENGINE CLASS ##
class Engine:
    # Upon declaration of a new engine:
    def __init__(self, thrust, P_inj, con_rat, MR=2, div_ang=15, con_ang=35, rad_rat=0.7, L_star=1.05):
        self.thrust = thrust # thrust: Design thrust [N]
        self.P_inj = P_inj # P_inj: Injector face pressure [Bar]
        self.con_rat = con_rat # con_rat: Enigne contraciton ratio
        self.MR = MR # Mixture ratio by weight (ox/fuel)
        self.div_ang = div_ang # Divergence half angle [Deg]
        self.con_ang = con_ang # Convergence half angle [Deg]
        self.rad_rat = rad_rat # rad_rat: R/Rmax (BETWEEN 0 & 1)
        self.L_star = L_star # Characteristic length [m]

    # Main design function
    def design_engine(self):
        # Call in global constants
        constants = Constants()

        # ==== Specify Parameters ====
        P_inj_psi = self.P_inj * constants.psi_to_bar # Psi is used in many CEA funcitons
        P_rat = self.P_inj/constants.P_amb # Inverse pressure ratio of overall expansion

        # Define RocketCEA object
        ispObj = CEA_Obj( oxName='LOX', fuelName='JetA', fac_CR=self.con_rat)

        # ==== Get Expansion Ratio ====

        # Throat pressure ratio (injector pressure / throat pressure)
        pip_t = ispObj.get_Throat_PcOvPe(P_inj_psi, self.MR)

        # Get the area ratio where optimal expansion happens
        self.expRat = ispObj.get_eps_at_PcOvPe(P_inj_psi, self.MR, P_rat)

        # ==== Get Throat Radius ====

        self.M = ispObj.get_MachNumber(P_inj_psi, self.MR, self.expRat) # Calculate mach number

        sonic_v = ispObj.get_SonicVelocities(P_inj_psi, self.MR, self.expRat) # Calculate sonic velocity

        self.sonic_v = sonic_v[2]*constants.ft_to_m # Convert from ft/s to m/s

        self.nozCorrFactor = (1 + math.cos(self.div_ang*math.pi/180))/2 # Correction factor for nozzle exit velocity, using conic appx
        self.V_exit = self.sonic_v * self.M # Sonic Vel * Mach @ Exit (m/s)

        self.mDot_tot = self.thrust/(self.nozCorrFactor * self.V_exit) # Calculate total mass flow based on V_exit (kg/s)
        self.mDot_o = (self.mDot_tot/(self.MR + 1)) * self.MR # Oxidizer mass flow (kg/s)
        self.mDot_f = self.mDot_tot - self.mDot_o # Fuel mass flow (kg/s)

        self.T = ispObj.get_Temperatures(P_inj_psi, self.MR, self.expRat, frozen=0, frozenAtThroat=0) # Get temperature array

        self.T_t = (self.T[1]-32)*5/9 + constants.C_to_K # Throat temperature (F to K)

        (self.MW, gamma) = ispObj.get_Throat_MolWt_gamma(P_inj_psi, self.MR) # Get molecular weight and specific heat ratio

        self.A_t = (self.mDot_tot/((self.P_inj / pip_t) * constants.bar)) * math.sqrt((constants.Ru * self.T_t)/(self.MW * gamma)) # Area at throat (m^2)

        self.R_t = math.sqrt(self.A_t / math.pi) # Calculate Throat Radius (m)

        self.con_ang = self.con_ang * math.pi / 180 # Nozzle contraction angle (degrees to rad)

        # Generate engine contour from external function:
        (self.engineContour, self.chBarrel, self.nozzleContour, self.R_tCurve, self.throatInd, self.conLeadInRadius) = getContour(self.R_t, self.L_star, self.con_rat, self.con_ang, self.div_ang, self.rad_rat, self.expRat)

        # Generate engine property array from helper function
        # This creates a 200x23 array describing various properties along the nozzle.
        # See getProperties.py for details
        self.engineProps = getProps(self.chBarrel, self.nozzleContour, self.throatInd, ispObj, P_inj_psi, self.MR, self.A_t)

        # Check for bad design (This occurs when thrust is way higher than chamber pressure should be and barrel becomes negative)
        if self.chBarrel[1,1] < 0:
            raise Exception("Thrust is too high for given chamber pressure to create a geometry with this method. \
                Try lowering thrust, increasing chamber pressure and/or decreasing contraction ratio.")

        # Performance parameters (not used in calculations)
        self.C_f = ispObj.get_PambCf(constants.psi_to_atm, P_inj_psi, self.MR, self.expRat)
        self.C_star = ispObj.get_Cstar(P_inj_psi, self.MR) * constants.ft_to_m * 0.975 # Calculate C*, m/s, apply 0.975 correction factor per H&H (pp. 70)


    # Compare CEA results with first order estimates from isentropic flow relations
    # Done mainly out of curiousity, but also validates results in case of mistakes
    def isentropic_comparison(self):

        # Throat Pressure Ratio
        gam = self.engineProps[0, 15] # Injector gamma
        pip_t_i = ((gam+1)/2)**(1/(gam-1))
        pip_t_cea = self.engineProps[99+self.throatInd, 2]
        print(self.throatInd)
        print(pip_t_i, pip_t_cea)

        # Exit Velocity
