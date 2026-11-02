import numpy as np
from rocketcea.cea_obj import CEA_Obj
import EngineParameters as ep


def sizeBaseEngine(D_c=None):
    """
    Computes 0D thermodynamic throat sizing, expansion ratio, and mass flow
    using RocketCEA and parameters from EngineParameters.py.
    """
    if D_c is None:
        D_c = ep.D_c

    rocket = CEA_Obj(oxName=ep.oxidizer, fuelName=ep.fuel)

    # --- Thermodynamic Sizing ---
    eps = rocket.get_eps_at_PcOvPe(ep.Pc, ep.MR, ep.Pc / ep.Pe, 1, 0)

    Isp_vac_ideal, Cstar_ideal, _ = rocket.get_IvacCstrTc(ep.Pc, ep.MR, eps, 1, 0)
    Cf_vac_ideal = (Isp_vac_ideal * 32.17405) / Cstar_ideal
    Cf_ideal = Cf_vac_ideal - (ep.Pe / ep.Pc) * eps

    Cstar_actual = Cstar_ideal * ep.Cstar_eff
    Cf_actual = Cf_ideal * ep.Cf_eff

    # Size the Throat Area (A_t) and Radius (R_t)
    A_t = ep.T_target / (Cf_actual * ep.Pc)  # in^2
    R_t = np.sqrt(A_t / np.pi)               # in

    # Chamber Contraction Ratio
    R_c = D_c / 2.0                          # in
    eps_c = (R_c / R_t)**2

    # Mass Flow Rate
    mdot_lbm = (ep.Pc * A_t * 32.17405) / Cstar_actual
    mdot_kg = mdot_lbm * 0.453592            # kg/s

    return {
        "A_t": A_t,
        "R_t": R_t,
        "R_c": R_c,
        "eps": eps,
        "eps_c": eps_c,
        "Cstar_actual": Cstar_actual,
        "Cf_actual": Cf_actual,
        "mdot_lbm": mdot_lbm,
        "mdot_kg": mdot_kg,
    }