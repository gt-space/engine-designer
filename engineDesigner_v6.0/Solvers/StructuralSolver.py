
import copy
import numpy as np
import EngineParameters as ep
from Geometry.LinerShell import getLinerGeometry, getOuterShellGeometry
from Functions.Structural import calcRevolvedShellMass
from Solvers.ChannelSolver import CoolantChannelSolver


def solveRegenStructuralShells(contour, channel, hot_gas, coolant_flow=None,
                               liner_material=None, shell_material=None,
                               couple_liner_heat_transfer=True, max_iter=2):
    """
    Coupled Thermo-Structural Solver:
    1. Evaluates Liner thickness (max of structural vs. thermal) via Geometry.LinerShell.
    2. Iterates with CoolantChannelSolver so variable t_liner(z) updates wall conduction
       resistance (R_wall = tw / (Ag * kw)) and channel pitch.
    3. Evaluates Outer Shell thickness (max of structural vs. thermal) after cooling channels.
    4. Integrates Liner, Rib, and Outer Shell masses using matprotlib densities.
    """
    channel_work = copy.deepcopy(channel)
    if liner_material is not None:
        channel_work["material"] = ep.load_matprotlib_material(liner_material)

    if coolant_flow is None:
        coolant_flow = CoolantChannelSolver(ep.T_cool_in, contour, channel_work, hot_gas)

    n_loops = max_iter if couple_liner_heat_transfer else 1
    for _ in range(n_loops):
        liner_geom = getLinerGeometry(
            contour, channel_work, hot_gas, coolant_flow,
            material=channel_work["material"]
        )
        tw_new = liner_geom["t_liner_selected"]

        if np.max(np.abs(tw_new - channel_work["tw"])) < 1e-5 or not couple_liner_heat_transfer:
            channel_work["tw"] = tw_new
            break

        # Update channel liner thickness and local pitch
        channel_work["tw"] = tw_new
        r_ch = contour["r"] + tw_new + 0.5 * channel_work["h_ch"]
        pitch = 2.0 * np.pi * r_ch / channel_work["N_ch"]
        channel_work["w_fin"] = np.maximum(pitch - channel_work["w_ch"], 0.2e-3)
        channel_work["f_fin"] = channel_work["w_fin"] / pitch

        coolant_flow = CoolantChannelSolver(ep.T_cool_in, contour, channel_work, hot_gas)

    shell_geom = getOuterShellGeometry(
        contour, channel_work, coolant_flow,
        t_liner_profile=channel_work["tw"],
        liner_material=channel_work["material"],
        shell_material=shell_material
    )

    rho_liner = float(liner_geom["material"].get("density"))
    rho_shell = float(shell_geom["material"].get("density"))

    _, mass_liner = calcRevolvedShellMass(
        contour["z"], contour["r"], channel_work["tw"], rho_liner, path_type="axial"
    )
    r_rib_in = contour["r"] + channel_work["tw"]
    _, mass_ribs_full = calcRevolvedShellMass(
        contour["z"], r_rib_in, channel_work["h_ch"], rho_liner, path_type="axial"
    )
    mass_ribs = mass_ribs_full * float(np.mean(channel_work["f_fin"]))

    _, mass_shell = calcRevolvedShellMass(
        contour["z"], shell_geom["r_shell_inner"], shell_geom["t_shell_selected"], rho_shell, path_type="axial"
    )

    return {
        "channel_updated": channel_work,
        "coolant_flow_updated": coolant_flow,
        "liner": liner_geom,
        "shell": shell_geom,
        "mass_liner_kg": mass_liner + mass_ribs,
        "mass_shell_kg": mass_shell,
        "mass_total_chamber_kg": mass_liner + mass_ribs + mass_shell
    }