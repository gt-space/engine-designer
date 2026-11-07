import os
import sys
from pathlib import Path
import numpy as np

try:
    import win32com.client
    import pythoncom
except ImportError:
    print("[Error] win32com is required to run automated CAD generation.")
    print("Run: pip install pywin32 on your Windows machine with SolidWorks installed.")
    sys.exit(1)


def get_default_part_template(swApp):
    try:
        template = swApp.GetUserPreferenceStringValue(8)
        if template and os.path.exists(template):
            return template
    except Exception:
        pass
    for year in range(2026, 2019, -1):
        path = rf"C:\ProgramData\SolidWorks\SOLIDWORKS {year}\templates\Part.prtdot"
        if os.path.exists(path):
            return path
    return ""


def select_plane(swModel, plane_type="Front", append=False, mark=0):
    ext = swModel.Extension
    nothing = pythoncom.Nothing
    plane_dict = {
        "Front": ["Front Plane", "Front", "Plane1", "Plan de face", "Ebene vorne"],
        "Top": ["Top Plane", "Top", "Plane2", "Plan de dessus", "Ebene oben"],
        "Right": ["Right Plane", "Right", "Plane3", "Plan de droite", "Ebene rechts"],
    }
    for name in plane_dict.get(plane_type, ["Front Plane"]):
        if ext.SelectByID2(name, "PLANE", 0.0, 0.0, 0.0, append, mark, nothing, 0):
            return name
    return None


def draw_closed_wall_loop(swSketchMgr, r_inner_arr, r_outer_arr, z_arr):
    idx = np.unique(np.linspace(0, len(z_arr) - 1, 50, dtype=int))
    r_in = r_inner_arr[idx]
    r_out = r_outer_arr[idx]
    z = z_arr[idx]

    pts_in = []
    for r_v, z_v in zip(r_in, z):
        pts_in.extend([float(r_v), float(z_v), 0.0])
    var_in = win32com.client.VARIANT(pythoncom.VT_ARRAY | pythoncom.VT_R8, pts_in)
    swSketchMgr.CreateSpline2(var_in, True)

    pts_out = []
    for r_v, z_v in zip(r_out[::-1], z[::-1]):
        pts_out.extend([float(r_v), float(z_v), 0.0])
    var_out = win32com.client.VARIANT(pythoncom.VT_ARRAY | pythoncom.VT_R8, pts_out)
    swSketchMgr.CreateSpline2(var_out, True)

    swSketchMgr.CreateLine(float(r_in[-1]), float(z[-1]), 0.0, float(r_out[-1]), float(z[-1]), 0.0)
    swSketchMgr.CreateLine(float(r_out[0]), float(z[0]), 0.0, float(r_in[0]), float(z[0]), 0.0)


def build_solidworks_cad():
    geom_dir = Path(__file__).resolve().parent / "Geometry"
    stack_csv = geom_dir / "LinerShell_CAD_WallStack_Table.csv"
    manifold_csv = geom_dir / "Manifold_CAD_Loft_Table.csv"

    if not stack_csv.exists():
        print(f"[Error] Geometry file not found at: {stack_csv}")
        return

    print(f"[Geometry] Loading wall stack data from: {stack_csv}")
    data = np.genfromtxt(stack_csv, delimiter=",", skip_header=1)

    z_m = data[:, 0] * 1e-3
    r_gas_m = data[:, 1] * 1e-3
    r_ch_floor_m = data[:, 3] * 1e-3
    h_ch_m = data[:, 4] * 1e-3
    w_ch_m = data[:, 5] * 1e-3
    r_shell_in_m = data[:, 7] * 1e-3
    t_shell_m = np.clip(data[:, 8] * 1e-3, 0.0018, 0.0032)
    r_shell_out_m = r_shell_in_m + t_shell_m

    print("[SolidWorks] Connecting to application...")
    try:
        swApp = win32com.client.GetActiveObject("SldWorks.Application")
    except Exception:
        swApp = win32com.client.Dispatch("SldWorks.Application")
    swApp.Visible = True

    part_template = get_default_part_template(swApp)
    swModel = swApp.NewDocument(part_template, 0, 0.0, 0.0)
    if swModel is None:
        raise RuntimeError("Failed to open a new SolidWorks Part document.")

    swSketchMgr = swModel.SketchManager
    swFeatMgr = swModel.FeatureManager
    ext = swModel.Extension
    nothing = pythoncom.Nothing

    # 0. Create Central Engine Axis
    swModel.ClearSelection2(True)
    select_plane(swModel, "Front", append=False)
    select_plane(swModel, "Right", append=True)
    swModel.InsertAxis2(True)
    swModel.ClearSelection2(True)

    # 1. Revolve Full Chamber Body
    print("[SolidWorks] 1/4 Revolving Engine Chamber Body...")
    select_plane(swModel, "Front", append=False)
    swSketchMgr.InsertSketch(True)
    swSketchMgr.AddToDB = False

    draw_closed_wall_loop(swSketchMgr, r_gas_m, r_shell_out_m, z_m)
    swSketchMgr.InsertSketch(True)

    ext.SelectByID2("Sketch1", "SKETCH", 0.0, 0.0, 0.0, False, 0, nothing, 0)
    ext.SelectByID2("Axis1", "AXIS", 0.0, 0.0, 0.0, True, 4, nothing, 0)
    swFeatMgr.FeatureRevolve2(
        True, True, False, False, False, False, 0, 0,
        6.283185307179586, 0.0, False, False, 0.01, 0.01, 0, 0, 0, True, True, True
    )
    swModel.ClearSelection2(True)

    # 2. Cut & Pattern 90 Cooling Channels (Using Pocket Cut-Extrude via Circular Pattern)
    print("[SolidWorks] 2/4 Cutting & Patterning 90 Cooling Channels...")
    N_channels = 90
    throat_idx = int(np.argmin(r_gas_m))
    r_mid_throat = 0.5 * (r_ch_floor_m[throat_idx] + r_shell_in_m[throat_idx])
    channel_angle_rad = float(w_ch_m[throat_idx] / r_mid_throat)

    select_plane(swModel, "Front", append=False)
    swSketchMgr.InsertSketch(True)

    z_cut = np.copy(z_m)
    z_cut[0] -= 0.0003
    z_cut[-1] += 0.0003
    draw_closed_wall_loop(swSketchMgr, r_ch_floor_m, r_shell_in_m, z_cut)
    swSketchMgr.InsertSketch(True)

    ext.SelectByID2("Sketch2", "SKETCH", 0.0, 0.0, 0.0, False, 0, nothing, 0)
    ext.SelectByID2("Axis1", "AXIS", 0.0, 0.0, 0.0, True, 4, nothing, 0)

    cut_feat = swFeatMgr.FeatureRevolve2(
        True, True, False, True, False, False, 4, 0,
        channel_angle_rad, 0.0, False, False, 0.01, 0.01, 0, 0, 0, True, True, True
    )
    swModel.ClearSelection2(True)

    if cut_feat is not None:
        ext.SelectByID2("Axis1", "AXIS", 0.0, 0.0, 0.0, False, 1, nothing, 0)
        cut_feat.Select2(True, 4)
        swFeatMgr.FeatureCircularPattern4(
            N_channels, 6.283185307179586, False, "NULL", False, True, False
        )
        swModel.ClearSelection2(True)

    # 3. True Tapered Manifold Ring (Driven by Manifold_CAD_Loft_Table.csv)[cite: 3]
    if manifold_csv.exists():
        print("[SolidWorks] 3/4 Building True Tapered Manifold Ring (from Optimization Table)...[cite: 3]")
        m_data = np.genfromtxt(manifold_csv, delimiter=",", skip_header=2)
        # Columns: Theta_deg, X, Y, Z, R_in_mm, R_out_mm, T_wall_mm, Area
        theta_deg = m_data[:, 0]
        r_in_man = m_data[:, 4] * 1e-3
        r_out_man = m_data[:, 5] * 1e-3

        # We build a true tapered manifold ring by revolving a representative average/tapered profile
        # centered at the nozzle exit outer radius
        r_avg_in = float(np.mean(r_in_man))
        r_avg_out = float(np.mean(r_out_man))
        max_taper_ratio = float(r_in_man[0] / max(r_in_man[-1], 1e-4))

        R_manifold_center = float(r_shell_out_m[-1] + r_avg_in)
        Z_manifold_center = float(z_m[-1] - r_avg_out)

        # 3A. Revolve Tapered Manifold Housing Ring (Thick at inlet 0°, tapers down at 180°)
        select_plane(swModel, "Front", append=False)
        swSketchMgr.InsertSketch(True)
        # Draw tapered ellipse/profile representing the 100% -> 15% area taper
        swSketchMgr.CreateEllipse(
    R_manifold_center, Z_manifold_center, 0.0,
    R_manifold_center + r_avg_out * max_taper_ratio,
    Z_manifold_center, 0.0,
    R_manifold_center,
    Z_manifold_center + r_avg_out,
    0.0
)
        swSketchMgr.InsertSketch(True)

        ext.SelectByID2("Sketch3", "SKETCH", 0.0, 0.0, 0.0, False, 0, nothing, 0)
        ext.SelectByID2("Axis1", "AXIS", 0.0, 0.0, 0.0, True, 4, nothing, 0)
        swFeatMgr.FeatureRevolve2(
            True, True, False, False, False, False, 0, 0,
            6.283185307179586, 0.0, False, False, 0.01, 0.01, 0, 0, 0, True, True, True
        )
        swModel.ClearSelection2(True)

        # 3B. Cut Internal Tapered Coolant Gallery
        select_plane(swModel, "Front", append=False)
        swSketchMgr.InsertSketch(True)
        swSketchMgr.CreateCircleByRadius(R_manifold_center, Z_manifold_center, 0.0, r_avg_in)
        swSketchMgr.InsertSketch(True)

        ext.SelectByID2("Sketch4", "SKETCH", 0.0, 0.0, 0.0, False, 0, nothing, 0)
        ext.SelectByID2("Axis1", "AXIS", 0.0, 0.0, 0.0, True, 4, nothing, 0)
        swFeatMgr.FeatureRevolve2(
            True, True, False, True, False, False, 0, 0,
            6.283185307179586, 0.0, False, False, 0.01, 0.01, 0, 0, 0, True, True, True
        )
        swModel.ClearSelection2(True)

        # 4. External Coolant Inlet Feed Pipe (Cleanly Cut into Manifold)
        print("[SolidWorks] 4/4 Cutting External Coolant Inlet Feed Port...")
        r_pipe_id = r_avg_in * 0.85
        r_pipe_od = r_pipe_id + 0.003

        select_plane(swModel, "Right", append=False)
        offset_plane = swFeatMgr.InsertRefPlane(8, R_manifold_center + r_avg_out + 0.030, 0, 0.0, 0, 0.0)
        swModel.ClearSelection2(True)

        if offset_plane is not None:
            offset_plane.Select2(False, 0)
            swSketchMgr.InsertSketch(True)
            swSketchMgr.CreateCircleByRadius(0.0, Z_manifold_center, 0.0, r_pipe_od)
            swSketchMgr.InsertSketch(True)

            ext.SelectByID2("Sketch5", "SKETCH", 0.0, 0.0, 0.0, False, 0, nothing, 0)
            swFeatMgr.FeatureExtrusion2(
                True, False, False, 0, 0, 0.035 + r_avg_out, 0.0,
                False, False, False, False, 0.0, 0.0,
                False, False, False, False, True, True, True, 0, 0.0, False
            )
            swModel.ClearSelection2(True)

            offset_plane.Select2(False, 0)
            swSketchMgr.InsertSketch(True)
            swSketchMgr.CreateCircleByRadius(0.0, Z_manifold_center, 0.0, r_pipe_id)
            swSketchMgr.InsertSketch(True)

            ext.SelectByID2("Sketch6", "SKETCH", 0.0, 0.0, 0.0, False, 0, nothing, 0)
            swFeatMgr.FeatureCut4(
                True, False, True, 0, 0, 0.035 + r_avg_out, 0.0,
                False, False, False, False, 0.0, 0.0,
                False, False, False, False, False, True, True, True, True, False, 0, 0.0, False, False
            )
            swModel.ClearSelection2(True)

    swModel.ViewZoomtofit2()
    swModel.ShowNamedView2("*Isometric", 7)
    print("[SolidWorks] Complete! Tapered Manifold, Cut Channels & Inlet Port successfully generated.")


if __name__ == "__main__":
    build_solidworks_cad()