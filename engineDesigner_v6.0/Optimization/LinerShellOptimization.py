import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize

from Solvers.StructuralSolver import solveRegenStructuralShells
from Functions.Structural import calcRevolvedShellMass

#Finds the minimum-mass smooth PCHIP spline thickness schedules for both the inner liner and outer closeout shell
def optimizeLinerAndShellProfiles(contour, channel, hot_gas, coolant_flow=None,
                                  liner_material=None, shell_material=None, n_ctrl=5):
    
    # 1. Solve baseline exact station-by-station requirements
    base_res = solveRegenStructuralShells(
        contour, channel, hot_gas, coolant_flow,
        liner_material=liner_material, shell_material=shell_material,
        couple_liner_heat_transfer=True
    )

    z = contour["z"]
    r_gas = contour["r"]
    r_shell_in = base_res["shell"]["r_shell_inner"]

    t_liner_req = base_res["liner"]["t_liner_selected"]
    t_liner_max = base_res["liner"]["t_liner_max_fourier"]
    t_shell_req = base_res["shell"]["t_shell_selected"]

    rho_liner = float(base_res["liner"]["material"].get("density"))
    rho_shell = float(base_res["shell"]["material"].get("density"))

    # Control stations anchored at chamber start, pre-throat, throat, post-throat, exit
    throat_idx = contour["throat_index"]
    ctrl_indices = np.unique(np.linspace(0, len(z) - 1, n_ctrl, dtype=int))
    if throat_idx not in ctrl_indices:
        ctrl_indices = np.sort(np.append(ctrl_indices[:-1], throat_idx))
    z_ctrl = z[ctrl_indices]

    def optimize_single_wall(t_req_array, r_base_array, rho_mat, t_upper_array=None):
        x0 = t_req_array[ctrl_indices] * 1.05
        bounds = [(float(t_req_array[idx]), float(t_upper_array[idx]) if t_upper_array is not None else 0.05)
                  for idx in ctrl_indices]

        def interp_wall(x_ctrl):
            return PchipInterpolator(z_ctrl, x_ctrl)(z)

        def objective(x_ctrl):
            t_prof = interp_wall(x_ctrl)
            _, mass = calcRevolvedShellMass(z, r_base_array, t_prof, rho_mat, path_type="axial")
            return mass

        constraints = [
            # Priority lower bound: t_smooth(z) >= max(t_struct(z), t_thermal(z)) everywhere
            {"type": "ineq", "fun": lambda x: (interp_wall(x) - t_req_array) * 1e3}
        ]
        if t_upper_array is not None:
            constraints.append(
                {"type": "ineq", "fun": lambda x: (t_upper_array - interp_wall(x)) * 1e3}
            )

        res = minimize(objective, x0, method="SLSQP", bounds=bounds, constraints=constraints)
        t_opt_prof = np.maximum(interp_wall(res.x), t_req_array)
        return t_opt_prof, res.x

    t_liner_smooth, liner_ctrl_pts = optimize_single_wall(
        t_liner_req, r_gas, rho_liner, t_upper_array=t_liner_max
    )
    t_shell_smooth, shell_ctrl_pts = optimize_single_wall(
        t_shell_req, r_shell_in, rho_shell
    )

    # Update dictionaries with smooth CAD-ready profiles
    base_res["liner"]["t_liner_selected"] = t_liner_smooth
    base_res["liner"]["r_liner_outer"] = r_gas + t_liner_smooth
    base_res["liner"]["z_ctrl"] = z_ctrl
    base_res["liner"]["t_ctrl_pts"] = liner_ctrl_pts

    base_res["shell"]["t_shell_selected"] = t_shell_smooth
    base_res["shell"]["r_shell_outer"] = r_shell_in + t_shell_smooth
    base_res["shell"]["z_ctrl"] = z_ctrl
    base_res["shell"]["t_ctrl_pts"] = shell_ctrl_pts

    base_res["channel_updated"]["tw"] = t_liner_smooth

    return base_res