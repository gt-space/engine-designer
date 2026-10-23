import numpy as np
from scipy.optimize import minimize

from Geometry import CoolingChannel
from Solvers import ChannelSolver

def setCoolantConstraints(T_cool_in, max_Twg, max_dP):
    return {
        "T_cool_in": T_cool_in,
        "Max_Twg": max_Twg,
        "Max_DP": max_dP
    }


def setGeometryConstraints(
    AR_bounds,
    w_min,
    w_max,
    fin_min,
    tw,
    N_ch,
    h_chamber_bounds,
    h_throat_bounds,
    h_exit_bounds,
    f_fin_chamber_bounds,
    f_fin_throat_bounds,
    f_fin_exit_bounds):

    return {
        "AR_Bounds": AR_bounds,
        "Min_W": w_min,
        "Max_W": w_max,
        "Min_Fin_W": fin_min,
        "tw": tw,
        "N_ch": N_ch,

        "h_chamber_bounds": h_chamber_bounds,
        "h_throat_bounds": h_throat_bounds,
        "h_exit_bounds": h_exit_bounds,

        "f_fin_chamber_bounds": f_fin_chamber_bounds,
        "f_fin_throat_bounds": f_fin_throat_bounds,
        "f_fin_exit_bounds": f_fin_exit_bounds}


def optimizeChannels(contour, hot_gas, Coolant_Constraints, Geometry_Constraints):

    # Coolant Constraints
    T_cool_in = Coolant_Constraints["T_cool_in"]
    max_Twg = Coolant_Constraints["Max_Twg"]
    max_DP = Coolant_Constraints["Max_DP"]

    #Geometry Constraints / Constants
    AR_min, AR_max = Geometry_Constraints["AR_Bounds"]
    min_W = Geometry_Constraints["Min_W"]
    max_W = Geometry_Constraints["Max_W"]
    min_fin_W = Geometry_Constraints["Min_Fin_W"]
    tw = Geometry_Constraints["tw"]
    N_ch = Geometry_Constraints["N_ch"]

    # Variable bounds, these are parameters that will be optimized
    bounds = [
        Geometry_Constraints["h_chamber_bounds"],
        Geometry_Constraints["h_throat_bounds"],
        Geometry_Constraints["h_exit_bounds"],

        Geometry_Constraints["f_fin_chamber_bounds"],
        Geometry_Constraints["f_fin_throat_bounds"],
        Geometry_Constraints["f_fin_exit_bounds"]]

    def buildChannel(x):

        h_chamber, h_throat, h_exit, f_fin_chamber, f_fin_throat, f_fin_exit = x

        channel = CoolingChannel.getCoolingChannelGeometry(
            contour,
            N_ch,
            tw,
            h_chamber,
            h_throat,
            h_exit,
            f_fin_chamber,
            f_fin_throat,
            f_fin_exit)

        # Reject impossible geometry
        if np.any(channel["w_ch"] <= 0):
            raise ValueError("Invalid geometry: channel width <= 0")
        if np.any(channel["w_fin"] <= 0):
            raise ValueError("Invalid geometry: fin width <= 0")
        if np.any(channel["h_ch"] <= 0):
            raise ValueError("Invalid geometry: channel height <= 0")
        if np.any(channel["A_ch"] <= 0):
            raise ValueError("Invalid geometry: channel area <= 0")
        if np.any(channel["Dh"] <= 0):
            raise ValueError("Invalid geometry: hydraulic diameter <= 0")

        return channel

    # --------------------------------
    # Solve one channel design
    # --------------------------------
    def evaluateDesign(x):

        try:
            channel = buildChannel(x)

            coolant_flow = ChannelSolver.CoolantChannelSolver(
                T_cool_in,
                contour,
                channel,
                hot_gas
            )

            # Reject NaNs / infinities
            for key in ["P_cool", "T_cool", "Twg", "Twc", "q_flux"]:
                if not np.all(np.isfinite(coolant_flow[key])):
                    raise ValueError(f"Invalid values in {key}")

            # Reject nonphysical pressures
            if np.any(coolant_flow["P_cool"] <= 0):
                raise ValueError("Coolant pressure became non-positive")

            return coolant_flow, True

        except (ValueError, RuntimeError, OverflowError):
            return  None, False


    # --------------------------------
    # Objective -> Minimize dP while maintaining temperature & geometry constraints
    # --------------------------------
    def objective(x):
        coolant_flow, valid = evaluateDesign(x)

        if not valid:
            return 1e6

        P_in = coolant_flow["P_cool"][-1]
        P_out = coolant_flow["P_cool"][0]

        dP = P_in - P_out

        return dP / 6894.75729


    # --------------------------------
    # Initial Guess
    # --------------------------------
    x0 = [
        6e-3,    # h_chamber
        2e-3,    # h_throat
        5e-3,    # h_exit

        0.40,    # fin fraction chamber
        0.40,    # fin fraction throat
        0.40     # fin fraction exit
    ]

     # --------------------------------
    # Constraint Functions
    # --------------------------------

    def constraint_Twg(x):
        coolant_flow, valid = evaluateDesign(x)

        if not valid:
            return -1e3
        
        return (max_Twg - np.max(coolant_flow["Twg"])) / max_Twg


    def constraint_dP(x):

        coolant_flow, valid = evaluateDesign(x)

        if not valid:
            return -1e6

        P_in = coolant_flow["P_cool"][-1]
        P_out = coolant_flow["P_cool"][0]

        dP = P_in - P_out

        return (max_DP - dP) / max_DP


    def constraint_minWidth(x):
        try:
            channel = buildChannel(x)
            return (np.min(channel["w_ch"]) - min_W) / min_W

        except ValueError:
            return -1e6
    def constraint_maxWidth(x):
        try:
            channel = buildChannel(x)

            return (max_W - np.max(channel["w_ch"])) / max_W

        except ValueError:
            return -1e3

    def constraint_minFinWidth(x):
        try:
            channel = buildChannel(x)
            return (np.min(channel["w_fin"]) - min_fin_W) / min_fin_W
        except ValueError:
            return -1e6
    
    def constraint_ARmin(x):
        try:
            channel = buildChannel(x)
            return (np.min(channel["AR_ch"]) - AR_min) / AR_min

        except ValueError:
            return -1e6


    def constraint_ARmax(x):
        try:
            channel = buildChannel(x)
            return (AR_max - np.max(channel["AR_ch"])) / AR_max

        except ValueError:
            return -1e6


    # --------------------------------
    # SciPy Constraints
    # --------------------------------
    scipy_constraints = [
        {"type": "ineq", "fun": constraint_Twg},
        {"type": "ineq", "fun": constraint_dP},
        {"type": "ineq", "fun": constraint_minWidth},
        {"type": "ineq", "fun": constraint_maxWidth},
        {"type": "ineq", "fun": constraint_minFinWidth},
        {"type": "ineq", "fun": constraint_ARmin},
        {"type": "ineq", "fun": constraint_ARmax}
    ]


    # --------------------------------
    # Run Optimization
    # --------------------------------
    result = minimize(
        objective,
        x0,
        method="SLSQP",
        bounds=bounds,
        constraints=scipy_constraints,
        options={
            "disp": True,
            "maxiter": 100})
    
    print("Optimization success:", result.success)
    print("Message:", result.message)
    print("Optimal variables:", result.x)


    # --------------------------------
    # Final Optimized Design
    # --------------------------------
    x_opt = result.x

    channel_optimized = buildChannel(x_opt)
    coolant_flow_optimized, valid = evaluateDesign(x_opt)

    if not valid:
        raise RuntimeError("Optimizer returned an invalid final design.")

    h_chamber,h_throat,h_exit,f_fin_chamber,f_fin_throat,f_fin_exit = x_opt


    # --------------------------------
    # Print optimized design variables
    # --------------------------------
    print(f"h_chamber = {h_chamber*1000:.3f} mm")
    print(f"h_throat  = {h_throat*1000:.3f} mm")
    print(f"h_exit    = {h_exit*1000:.3f} mm")

    print(f"f_fin_chamber = {f_fin_chamber:.3f}")
    print(f"f_fin_throat  = {f_fin_throat:.3f}")
    print(f"f_fin_exit    = {f_fin_exit:.3f}")


    # --------------------------------
    # Print resulting geometry
    # --------------------------------
    print(f"Min Channel Width = "
        f"{np.min(channel_optimized['w_ch'])*1000:.3f} mm")

    print(f"Min Fin Width = "
        f"{np.min(channel_optimized['w_fin'])*1000:.3f} mm")

    print(f"Max Fin Width = "
        f"{np.max(channel_optimized['w_fin'])*1000:.3f} mm")

    print(f"AR Range = "
        f"{np.min(channel_optimized['AR_ch']):.3f} - "
        f"{np.max(channel_optimized['AR_ch']):.3f}")


    # --------------------------------
    # Print performance
    # --------------------------------
    print(f"Max Twg = {np.max(coolant_flow_optimized['Twg']):.1f} K")
    dP_opt = (coolant_flow_optimized["P_cool"][-1]- coolant_flow_optimized["P_cool"][0])

    print(
        f"Pressure Drop = "
        f"{dP_opt / 6894.75729:.2f} psi")

    return {
        "result": result,
        "x_opt": x_opt,
        "channel": channel_optimized,
        "coolant_flow": coolant_flow_optimized
    }