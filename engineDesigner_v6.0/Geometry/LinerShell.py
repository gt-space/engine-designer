import numpy as np
from pathlib import Path
import EngineParameters as ep
from EngineParameters import load_matprotlib_material
from Functions.LinerShell import (
    calcLinerStructuralThickness,
    calcLinerThermalThickness,
    calcShellStructuralThickness,
    calcShellThermalThickness,
)
from Functions.Structural import calcRevolvedShellMass

#Evaluates structural and thermal liner thicknesses along the axial contour and selects the HIGHER thickness at every station.
def getLinerGeometry(contour, channel, hot_gas, coolant_flow,
                     material=None, safety_factor_struct=None,
                     safety_factor_thermal=None, t_min_mfg=None):
    
    mat = load_matprotlib_material(material if material is not None else channel.get("material", ep.liner_material_name))
    sf_s = safety_factor_struct if safety_factor_struct is not None else ep.SF_structural
    sf_t = safety_factor_thermal if safety_factor_thermal is not None else ep.SF_thermal
    t_min = t_min_mfg if t_min_mfg is not None else ep.t_liner_min_mfg

    r_gas = contour["r"]
    w_ch = channel["w_ch"]
    P_gas = hot_gas["P"]
    P_cool = coolant_flow["P_cool"]
    Twg = coolant_flow["Twg"]
    Twc = coolant_flow["Twc"]
    q_flux = coolant_flow["q_flux"]
    Tw_avg = 0.5 * (Twg + Twc)

    t_struct = calcLinerStructuralThickness(
        P_cool, P_gas, r_gas, w_ch, mat, Tw_avg,
        safety_factor=sf_s, t_min_mfg=t_min
    )
    t_thermal, t_max_fourier = calcLinerThermalThickness(
        P_cool, P_gas, w_ch, q_flux, Twg, Twc, mat,
        safety_factor=sf_t, t_min_mfg=t_min
    )

    # Higher thickness takes priority
    t_selected = np.maximum(t_struct, t_thermal)
    governing_mode = np.where(t_thermal >= t_struct, "Thermal", "Structural")

    return {
        "r_liner_inner": r_gas,
        "r_liner_outer": r_gas + t_selected,
        "t_liner_selected": t_selected,
        "t_liner_structural": t_struct,
        "t_liner_thermal": t_thermal,
        "t_liner_max_fourier": t_max_fourier,
        "governing_mode": governing_mode,
        "exceeds_fourier_limit": t_selected > t_max_fourier,
        "material": mat,
        "material_name": mat.name
    }

#Evaluates structural and thermal outer closeout shell thicknesses and selects the HIGHER thickness at every station.
def getOuterShellGeometry(contour, channel, coolant_flow, t_liner_profile=None,
                          liner_material=None, shell_material=None,
                          safety_factor_struct=None, safety_factor_thermal=None,
                          t_min_mfg=None):

    mat_liner = load_matprotlib_material(liner_material if liner_material is not None else channel.get("material", ep.liner_material_name))
    mat_shell = load_matprotlib_material(shell_material if shell_material is not None else ep.shell_material_name)
    sf_s = safety_factor_struct if safety_factor_struct is not None else ep.SF_structural
    sf_t = safety_factor_thermal if safety_factor_thermal is not None else ep.SF_thermal
    t_min = t_min_mfg if t_min_mfg is not None else ep.t_shell_min_mfg

    tw = t_liner_profile if t_liner_profile is not None else channel["tw"]
    r_shell_inner = contour["r"] + tw + channel["h_ch"]

    P_cool = coolant_flow["P_cool"]
    T_cool = coolant_flow["T_cool"]
    Twg = coolant_flow["Twg"]
    Twc = coolant_flow["Twc"]

    t_struct = calcShellStructuralThickness(
        P_cool, r_shell_inner, channel["w_ch"], mat_shell, T_cool,
        safety_factor=sf_s, t_min_mfg=t_min
    )
    t_thermal = calcShellThermalThickness(
        P_cool, r_shell_inner, tw, channel["h_ch"], channel["f_fin"],
        Twg, Twc, T_cool, mat_liner, mat_shell,
        safety_factor=sf_t, t_min_mfg=t_min
    )

    # Higher thickness takes priority
    t_selected = np.maximum(t_struct, t_thermal)
    governing_mode = np.where(t_thermal >= t_struct, "Thermal", "Structural")

    return {
        "r_shell_inner": r_shell_inner,
        "r_shell_outer": r_shell_inner + t_selected,
        "t_shell_selected": t_selected,
        "t_shell_structural": t_struct,
        "t_shell_thermal": t_thermal,
        "governing_mode": governing_mode,
        "material": mat_shell,
        "shell_material_name": mat_shell.name
    }

#Exports axial and radial coordinates (in mm) of the hot-gas contour, liner outer wall, channel roof / shell inner wall, and outer closeout shell wall
def exportLinerShellCADToCSV(contour, channel, liner_geom, shell_geom,
                             filename=None):
    if filename is None:
        geometry_dir = Path(__file__).resolve().parent
        output_path = geometry_dir / "LinerShell_CAD_WallStack_Table.csv"
    else:
        output_path = Path(filename)
        if not output_path.is_absolute():
            # If a relative path was passed, resolve it relative to Geometry/
            output_path = Path(__file__).resolve().parent / output_path.name

    # Ensure the parent folder exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    z_mm = contour["z"] * 1e3
    r_gas_mm = contour["r"] * 1e3
    t_liner_mm = liner_geom["t_liner_selected"] * 1e3
    r_liner_out_mm = liner_geom["r_liner_outer"] * 1e3
    r_shell_in_mm = shell_geom["r_shell_inner"] * 1e3
    t_shell_mm = shell_geom["t_shell_selected"] * 1e3
    r_shell_out_mm = shell_geom["r_shell_outer"] * 1e3
    w_ch_mm = channel["w_ch"] * 1e3
    h_ch_mm = channel["h_ch"] * 1e3
    w_fin_mm = channel["w_fin"] * 1e3

    header = (
        f"Regen Engine 3D CAD Wall Stack-Up Table | "
        f"Liner: {liner_geom['material_name']} | Shell: {shell_geom['shell_material_name']}\n"
        "Z_Axial_mm,R_HotGas_Inner_mm,T_Liner_mm,R_Channel_Floor_mm,"
        "Channel_Height_mm,Channel_Width_mm,Fin_Width_mm,"
        "R_Shell_Inner_mm,T_OuterShell_mm,R_Shell_Outer_mm"
    )

    table = np.column_stack((
        z_mm, r_gas_mm, t_liner_mm, r_liner_out_mm,
        h_ch_mm, w_ch_mm, w_fin_mm,
        r_shell_in_mm, t_shell_mm, r_shell_out_mm
    ))

    np.savetxt(str(output_path), table, delimiter=",", header=header, comments="# ", fmt="%.5f")
    print(f"[Geometry] Saved Liner & Shell CAD Wall Stack-Up Table to: {output_path}")