# DESIGN.PY – Defines engine class and applies relevant equations to build
# and simulate an engine geometry given input parameters

# IMPORTANT RESULTS:
# engineContour is a numpy array with the following column format:
# [R, Z]
# engineProps is a numpy array with the following column format:
# [R, Z, pip, aeat, mach, cf, ivac, isp, p_OW, t, rho_OW, h, u, mw, cp, gam, son,
# vis, cond, pran, cpfz, condfz, pranfz]

# engineProps Explanation:
# engine.engineProps[axial station number (0 to 200 by default), engine prop index]
# E.g engine.engineProps[0, 0] is initial CC radius
# E.g engine.engineProps[-1, 0] is initial CC radius

# Engine Properties Index:

#0 R
#1 Z
#2 pip
#3 aeat
#4 mach
#5 cf
#6 ivac
#7 isp
#8 p_OW
#9 t
#10 rho_OW
#11 h
#12 u
#13 mw
#14 cp
#15 gam
#16 son
#17 vis
#18 cond
#19 pran
#20 cpfz
#21 condfz
#22 pranfz (or -1)


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
    def __init__(self, thrust, P_inj, P_e, con_rat, adv_data, MR=2, wall_MR = 2, L_star=1.05, cstar_eff = 1, numPTS = 200, fuel="JetA",ox="LOX"):
        self.noz_correction = []
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
        self.ox = ox
        self.fuel = fuel

    # Main design function
    def design_engine(self):
        # Call in global constants
        constants = Constants()

        # Specify Parameters
        self.P_inj_psi = self.P_inj * constants.psi_to_bar # Psi is used in many CEA funcitons
        P_rat_tot = self.P_inj / self.P_e # Inverse pressure ratio of overall expansion

        # Define RocketCEA object
        con_rat=self.con_rat
        if self.con_rat==[]:
            con_rat= 4 #default contraction ratio value; use on the first iteration

        ispObj = CEA_Obj(oxName=self.ox, fuelName=self.fuel, fac_CR=con_rat)

        # Get Expansion Ratio
        pip_t = ispObj.get_Throat_PcOvPe(self.P_inj_psi, self.MR) # Throat pressure ratio (inj pressure / throat pressure)
        exp_rat = ispObj.get_eps_at_PcOvPe(self.P_inj_psi, self.MR, P_rat_tot)

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
        self.noz_correction = noz_correction

        # Get more params needed for throat sizing
        temps = ispObj.get_Temperatures(self.P_inj_psi, self.MR, exp_rat, frozen=0, frozenAtThroat=0) # Get temperature array
        T_t = (temps[1] - 32) * 5/9 + constants.C_to_K # Throat temperature [F to K]
        (MW_t, gamma) = ispObj.get_Throat_MolWt_gamma(self.P_inj_psi, self.MR) # Get molecular weight and specific heat ratio

        # Find throat area givens this equation: http://www.braeunig.us/space/propuls.htm [See 1.26]
        self.A_t = (self.mDot_tot/((self.P_inj / pip_t) * constants.bar)) * math.sqrt((constants.Ru * T_t)/(MW_t * gamma)) # Area at throat (m^2)

        self.R_t = math.sqrt(self.A_t / math.pi) # Calculate Throat Radius [m]
        if self.con_rat==[]:
            D_t_cm = self.R_t * 200
            self.con_rat = 8.0 * D_t_cm ** (-0.6) + 1.25
        # Generate engine contour from external function:
        (self.engineContour, chBarrel, nozzleContour, self.R_tCurve, self.throatInd, self.conLeadInRadius, self.theta_e, self.throat_end_ind, self.throatInd_engprops,self.chBarrel_endInd) = get_contour(self.R_t, self.con_rat, exp_rat, self.L_star, self.adv_data, self.numPTS)
        self.nozzleContour = nozzleContour
        # Generate engine property array from helper function
        # This creates a 200x23 array describing various properties along the nozzle.
        # See CEA_properties.py for details
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
        pip_t_cea = self.engineProps[self.numPTS + self.throatInd, 2]
        print("Isentropic Critical Pressure Ratio: " + str(pip_t_i))
        print("CEA Critical Pressure Ratio: " + str(pip_t_cea))
        print(" ~~~~~~~~~~~~~~~~~~~ ")

        # Exit Velocity

        # Mass Flow Rate
