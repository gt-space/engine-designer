# DESIGN.PY – Defines engine class and applies relevant equations to build
# and simulate an engine geometry given input parameters

# IMPORTANT RESULTS:
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

## DEFINE ENGINE CLASS ##
class Engine:
    # Upon declaration of a new engine:
    def __init__(self, thrust, P_inj, P_e, con_rat, adv_data, MR=2, wall_MR = 2, L_star=1.05, cstar_eff = 1, numPTS = 200):
        self.thrust = thrust # thrust: Design thrust [N]
        self.P_inj = P_inj # P_inj: Injector face pressure [bar]
        self.P_e = P_e # P_e: Exit pressure [bar]
        self.con_rat = con_rat # con_rat: Enigne contraciton ratio
        self.adv_data = adv_data # Geometric params for bell nozzle {dict}
        self.MR = MR # Core mixture ratio by weight (ox/fuel)
        self.wall_MR = wall_MR # Mixture ratio near wall by weight (ox/fuel)
        self.L_star = L_star # Characteristic length [m]
        self.cstar_eff = cstar_eff #C* efficiency
        self.numPTS = numPTS

    # Main design function
    def design_engine(self):
        # Call in global constants
        constants = Constants()

        # Specify Parameters
        self.P_inj_psi = self.P_inj * constants.psi_to_bar # Psi is used in many CEA funcitons
        P_rat_tot = self.P_inj / self.P_e # Inverse pressure ratio of overall expansion

        # Define RocketCEA object
        ispObj = CEA_Obj( oxName='LOX', fuelName='JetA', fac_CR=self.con_rat)

        # Get Expansion Ratio
        pip_t = ispObj.get_Throat_PcOvPe(self.P_inj_psi, self.MR) # Throat pressure ratio (inj pressure / throat pressure)
        exp_rat = ispObj.get_eps_at_PcOvPe(self.P_inj_psi, self.MR, P_rat_tot) # Get area ratio for desired expansion P_rat

        # Get Throat Radius
        # First need exit velocity
        mach_e = ispObj.get_MachNumber(self.P_inj_psi, self.MR, exp_rat) # Calculate exit mach number
        v_son_vec = ispObj.get_SonicVelocities(self.P_inj_psi, self.MR, exp_rat) # Calculate sonic velocities
        v_son_e = v_son_vec[2] * constants.ft_to_m # Get exit sonic velocity value @ exit. Convert from ft/s to m/s
        noz_correction = (1 + math.cos(self.adv_data["div_ang"] * math.pi/180))/2 # Correction factor for exit velocity, conical
        self.V_exit = v_son_e * mach_e # Exit velocity [m/s]

        # Then find mass flow rate given design thrust. Ignore pressure difference compofrom exhaust and ambient
        self.mDot_tot = self.thrust/(noz_correction * self.V_exit) # Calculate total mass flow based on V_exit [kg/s]
        self.mDot_o = (self.mDot_tot/(self.MR + 1)) * self.MR # Oxidizer mass flow [kg/s]
        self.mDot_f = self.mDot_tot - self.mDot_o # Fuel mass flow [kg/s]

        # Get more params needed for throat sizing
        temps = ispObj.get_Temperatures(self.P_inj_psi, self.MR, exp_rat, frozen=0, frozenAtThroat=0) # Get temperature array
        T_t = (temps[1] - 32) * 5/9 + constants.C_to_K # Throat temperature [F to K]
        (MW_t, gamma) = ispObj.get_Throat_MolWt_gamma(self.P_inj_psi, self.MR) # Get molecular weight and specific heat ratio

        # Find throat area given this equation: http://www.braeunig.us/space/propuls.htm [See 1.26]
        self.A_t = (self.mDot_tot/((self.P_inj / pip_t) * constants.bar)) * math.sqrt((constants.Ru * T_t)/(MW_t * gamma)) # Area at throat (m^2)
        self.R_t = math.sqrt(self.A_t / math.pi) # Calculate Throat Radius [m]

        # Generate engine contour from external function:
        (self.engineContour, chBarrel, nozzleContour, self.R_tCurve, self.throatInd, self.conLeadInRadius, self.theta_e, self.throat_end_ind) = get_contour(self.R_t, self.con_rat, exp_rat, self.L_star, self.adv_data, self.numPTS)

        # Generate engine property array from helper function
        # This creates a 200x23 array describing various properties along the nozzle.
        # See getProperties.py for details
        self.engineProps = get_props(chBarrel, nozzleContour, self.throatInd, ispObj, self.P_inj_psi, self.MR, self.A_t)
        self.wallProps = get_props(chBarrel, nozzleContour, self.throatInd, ispObj, self.P_inj_psi, self.wall_MR, self.A_t)

        # Check for bad design (This occurs when thrust is way higher than chamber pressure should be and barrel becomes negative)
        if chBarrel[1,1] < 0:
            raise Exception("Thrust is too high for given chamber pressure to create a geometry with this method. \
                Try lowering thrust, increasing chamber pressure and/or decreasing contraction ratio.")

        # Performance parameters (not used in calculations)
        self.C_f = ispObj.get_PambCf(constants.psi_to_atm, self.P_inj_psi, self.MR, exp_rat)
        self.C_star = ispObj.get_Cstar(self.P_inj_psi, self.MR) * constants.ft_to_m * self.cstar_eff# Calculate C*, m/s, apply 0.975 correction factor per H&H (pp. 70)
        self.mDot_f = self.mDot_f / self.cstar_eff


    # Compare CEA results with first order estimates from isentropic flow relations
    # Done mainly out of curiousity, but also validates results in case of mistakes
    def isentropic_comparison(self):
        print("*** ISENTROPIC COMPARISON TO CEA ***")
        # Throat Pressure Ratio
        gam = self.engineProps[0, 15] # Injector gamma
        pip_t_i = ((gam+1)/2)**(1/(gam-1)) # Critical pressure ratio [P_c/P_t]
        pip_t_cea = self.engineProps[100 + self.throatInd, 2]
        print("Isentropic Critical Pressure Ratio: " + str(pip_t_i))
        print("CEA Critical Pressure Ratio: " + str(pip_t_cea))
        print(" ~~~~~~~~~~~~~~~~~~~ ")

        # Exit Velocity

        # Mass Flow Rate
