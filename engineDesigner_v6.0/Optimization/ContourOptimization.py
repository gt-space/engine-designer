import numpy as np
from scipy.optimize import minimize
from rocketisp.geometry import Geometry

import EngineParameters as ep
from Functions.EngineSizing import sizeBaseEngine
from Geometry.Contour import search_Lc, evaluateContourGeometry, saveOptimizedContour


def setContourWeights(
    surface_area=5.0,
    engine_weight=20.0,
    price=10.0,
    line_loss=0.3,
    volume=0.9,
    manufacturability=50.0,
    thrust_penalty=100.0,
    cd_penalty=10.0
):
    # Returns a dictionary of cost function weights for easy tuning.
    return {
        "a_surf": surface_area,
        "b_weight": engine_weight,
        "c_price": price,
        "d_loss": line_loss,
        "e_vol": volume,
        "f_mfg": manufacturability,
        "g_thrust": thrust_penalty,
        "h_cd": cd_penalty,
    }


def setContourBounds(
    theta_c_range=(30.0, 45.0),
    Lstar_range=(40.0, 50.0),
    Lf_range=(60.0, 100.0),
    Rtc_range=(0.5, 1.5),
    D_c=6.5,
    Rtd_Rt=0.382
):
    # Returns geometric bounds and fixed contour constants.
    return {
        "theta_c_range": theta_c_range,
        "Lstar_range": np.array(Lstar_range),
        "Lf_range": Lf_range,
        "Rtc_range": Rtc_range,
        "D_c": D_c,
        "Rtd_Rt": Rtd_Rt,
    }




def optimizeContour(weights=None, geo_bounds=None):
    # Runs SLSQP contour optimization and exports the resulting contour.
    if weights is None:
        weights = setContourWeights()
    if geo_bounds is None:
        geo_bounds = setContourBounds()

    # 1. Size throat and expansion ratio
    sizing = sizeBaseEngine(geo_bounds["D_c"])
    A_t = sizing["A_t"]
    R_t = sizing["R_t"]
    eps = sizing["eps"]
    eps_c = sizing["eps_c"]
    Rtd_Rt = geo_bounds["Rtd_Rt"]

    # 2. Compute Lc bounds from L* bounds
    Vc_range = A_t * geo_bounds["Lstar_range"]
    Rtc_range = geo_bounds["Rtc_range"]
    theta_c_range = geo_bounds["theta_c_range"]
    Lf_range = geo_bounds["Lf_range"]

    G1 = Geometry(Rthrt=R_t, CR=eps_c, RupThroat=Rtc_range[1],
                  RchmConv=Rtc_range[1], cham_conv_deg=theta_c_range[0], LchamberInp=10)
    G2 = Geometry(Rthrt=R_t, CR=eps_c, RupThroat=Rtc_range[0],
                  RchmConv=Rtc_range[0], cham_conv_deg=theta_c_range[1], LchamberInp=10)

    Lc_range = np.around([search_Lc(Vc_range[0], G1), search_Lc(Vc_range[1], G2)], 1)

    bounds = [
        (Lc_range[0], Lc_range[1]),
        theta_c_range,
        Lf_range,
        Rtc_range
    ]

    # 3. Define Cost Function using the weights 
    def cost_function(x):
        Lc, theta_c_val, Lf_val, Rtc_val = x
        vol, surf = evaluateContourGeometry(R_t, eps, eps_c, Rtd_Rt, Lc, theta_c_val, Lf_val, Rtc_val)

        weight_pen = 5.0 * Lc if Lc > 10.0 else 0.0
        manufacturability = 1.0
        price = 2.0 * Lc
        line_loss = 1.0 + 0.54 / (eps_c**2)
        thrust_penalty = (100.0 - Lf_val) / 100.0
        cd_penalty = (3.0 - Rtc_val) / 3.0

        numerator = (
            weights["a_surf"] * surf
            + weights["b_weight"] * weight_pen
            + weights["c_price"] * price
            + weights["d_loss"] * line_loss
            + weights["g_thrust"] * thrust_penalty
            + weights["h_cd"] * cd_penalty
        )
        denominator = weights["e_vol"] * vol + weights["f_mfg"] * manufacturability
        return numerator / denominator

    # 4. Run Optimizer
    x0 = np.array([8.0, 37.5, 80.0, 1.5])
    result = minimize(cost_function, x0, bounds=bounds, method="SLSQP")
    Lc_opt, theta_c_opt, Lf_opt, Rtc_opt = np.around(result.x, 1)

    # 5. Save and Return
    params = saveOptimizedContour(R_t, eps, eps_c, Rtd_Rt, Lc_opt, theta_c_opt, Lf_opt, Rtc_opt)
    return params