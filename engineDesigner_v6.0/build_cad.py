import os
import re
import sys
from pathlib import Path
import numpy as np

try:
    import win32com.client
    import pythoncom
except ImportError:
    print("[Error] pywin32 is required to run automated CAD generation.")
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


def draw_closed_wall_loop(swSketchMgr, r_inner_arr, r_outer_arr, z_arr, n_pts=60):
    """Draws a watertight closed spline loop on the active sketch plane."""
    idx = np.unique(np.linspace(0, len(z_arr) - 1, n_pts, dtype=int))
    r_in = r_inner_arr[idx]
    r_out = r_outer_arr[idx]
    z = z_arr[idx]

    # Direct database insertion prevents SolidWorks auto-snapping bugs
    swSketchMgr.AddToDB = True

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

    # Close top (nozzle exit) and bottom (injector end) caps
    swSketchMgr.CreateLine(
        float(r_in[-1]), float(z[-1]), 0.0,
        float(r_out[-1]), float(z[-1]), 0.0
    )
    swSketchMgr.CreateLine(
        float(r_out[0]), float(z[0]), 0.0,
        float(r_in[0]), float(z[0]), 0.0
    )
    swSketchMgr.AddToDB = False


def color_latest_body(swModel, rgb_tuple, transparency=0.0):
    """Applies RGB color and optional transparency to the most recently created solid body."""
    try:
        bodies = swModel.GetBodies2(0, False)  # 0 = swSolidBody
        if bodies and len(bodies) > 0:
            body = bodies[-1]
            # [R, G, B, Ambient, Diffuse, Specular, Shininess, Transparency, Emission]
            mat_props = [
                rgb_tuple[0], rgb_tuple[1], rgb_tuple[2],
                0.5, 0.8, 0.5, 0.4, transparency, 0.0
            ]
            body.MaterialPropertyValues2 = win32com.client.VARIANT(
                pythoncom.VT_ARRAY | pythoncom.VT_R8, mat_props
            )
    except Exception:
        pass


def parse_inlet_angle_from_csv(manifold_csv_path):
    """Extracts the optimized inlet wye half-angle from the CSV header comment."""
    try:
        with open(manifold_csv_path, "r") as f:
            first_line = f.readline()
        match = re.search(r"Half-Angle:\s*([0-9.]+)", first_line)
        if match:
            return float(match.group(1))
    except Exception:
        pass
    return 45.0


def build_solidworks_cad():
    geom_dir = Path(__file__).resolve().parent / "Geometry"
    stack_csv = geom_dir / "LinerShell_CAD_WallStack_Table.csv"
    manifold_csv = geom_dir / "Manifold_CAD_Loft_Table.csv"

    if not stack_csv.exists():
        print(f"[Error] Geometry file not found at: {stack_csv}")
        print("Run `python main.py` first to generate the CAD tables.")
        return

    print(f"[Geometry] Loading wall stack data from: {stack_csv}")
    data = np.genfromtxt(stack_csv, delimiter=",", comments="#")

    z_m = data[:, 0] * 1e-3
    r_gas_m = data[:, 1] * 1e-3
    r_ch_floor_m = data[:, 3] * 1e-3
    w_ch_m = data[:, 5] * 1e-3
    r_shell_in_m = data[:, 7] * 1e-3
    t_shell_m = np.clip(data[:, 8] * 1e-3, 0.0015, 0.0060)
    r_shell_out_m = r_shell_in_m + t_shell_m

    print("[SolidWorks] Connecting to SolidWorks application...")
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

    # -------------------------------------------------------------------------
    # 0. Create Central Engine Axis (Axis1 = Intersection of Front & Right)
    # -------------------------------------------------------------------------
    swModel.ClearSelection2(True)
    select_plane(swModel, "Front", append=False)
    select_plane(swModel, "Right", append=True)
    swModel.InsertAxis2(True)
    swModel.ClearSelection2(True)

    # -------------------------------------------------------------------------
    # 1. Revolve Inner Liner + Fin Block (r_gas -> r_shell_in) [Copper Body]
    # -------------------------------------------------------------------------
    print("[SolidWorks] 1/5 Revolving Inner Hot-Gas Liner & Rib Block (GrCop-42)...")
    select_plane(swModel, "Front", append=False)
    swSketchMgr.InsertSketch(True)
    draw_closed_wall_loop(swSketchMgr, r_gas_m, r_shell_in_m, z_m)
    liner_sketch_name = swModel.GetActiveSketch2().Name
    swSketchMgr.InsertSketch(True)

    ext.SelectByID2(liner_sketch_name, "SKETCH", 0.0, 0.0, 0.0, False, 0, nothing, 0)
    ext.SelectByID2("Axis1", "AXIS", 0.0, 0.0, 0.0, True, 16, nothing, 0)
    swFeatMgr.FeatureRevolve2(
        True, True, False, False, False, False, 0, 0,
        2.0 * np.pi, 0.0, False, False, 0.01, 0.01, 0, 0, 0, True, True, True
    )
    swModel.ClearSelection2(True)
    color_latest_body(swModel, (0.85, 0.45, 0.20), transparency=0.0)  # Copper color

    # -------------------------------------------------------------------------
    # 2. Cut & Pattern 90 Cooling Channels into the Liner Body
    # -------------------------------------------------------------------------
    print("[SolidWorks] 2/5 Cutting & Patterning 90 Cooling Channels...")
    N_channels = 90
    throat_idx = int(np.argmin(r_gas_m))
    r_mid_throat = 0.5 * (r_ch_floor_m[throat_idx] + r_shell_in_m[throat_idx])
    channel_angle_rad = float(w_ch_m[throat_idx] / r_mid_throat)

    select_plane(swModel, "Front", append=False)
    swSketchMgr.InsertSketch(True)

    # Extend axial ends by 0.5 mm and radial roof by +0.2 mm so the open top of
    # the channel cleanly breaks out of the Liner body before the Shell is added
    z_cut = np.copy(z_m)
    z_cut[0] -= 0.0005
    z_cut[-1] += 0.0005
    r_cut_outer = r_shell_in_m + 0.0002

    draw_closed_wall_loop(swSketchMgr, r_ch_floor_m, r_cut_outer, z_cut)
    cut_sketch_name = swModel.GetActiveSketch2().Name
    swSketchMgr.InsertSketch(True)

    ext.SelectByID2(cut_sketch_name, "SKETCH", 0.0, 0.0, 0.0, False, 0, nothing, 0)
    ext.SelectByID2("Axis1", "AXIS", 0.0, 0.0, 0.0, True, 16, nothing, 0)

    # Dir1Type = 1 (swEndCondMidPlane) revolves symmetrically around the Front Plane
    cut_feat = swFeatMgr.FeatureRevolve2(
        True, True, False, True, False, False, 1, 0,
        channel_angle_rad, 0.0, False, False, 0.01, 0.01, 0, 0, 0, True, True, True
    )
    swModel.ClearSelection2(True)

    if cut_feat is not None:
        ext.SelectByID2("Axis1", "AXIS", 0.0, 0.0, 0.0, False, 1, nothing, 0)
        cut_feat.Select2(True, 4)
        # GeometryPattern = True speeds up 90-channel patterning significantly
        swFeatMgr.FeatureCircularPattern4(
            N_channels, 2.0 * np.pi, False, "NULL", True, True, False
        )
        swModel.ClearSelection2(True)
    else:
        print("[Warning] Single cooling channel revolve cut failed to create.")

    # -------------------------------------------------------------------------
    # 3. Revolve Outer Closeout Shell as a Separate Body (Merge = False)
    # -------------------------------------------------------------------------
    print("[SolidWorks] 3/5 Revolving Outer Closeout Shell (Separate Multi-Body)...")
    select_plane(swModel, "Front", append=False)
    swSketchMgr.InsertSketch(True)
    draw_closed_wall_loop(swSketchMgr, r_shell_in_m, r_shell_out_m, z_m)
    shell_sketch_name = swModel.GetActiveSketch2().Name
    swSketchMgr.InsertSketch(True)

    ext.SelectByID2(shell_sketch_name, "SKETCH", 0.0, 0.0, 0.0, False, 0, nothing, 0)
    ext.SelectByID2("Axis1", "AXIS", 0.0, 0.0, 0.0, True, 16, nothing, 0)
    # Merge = False (3rd boolean from end) keeps Shell distinct from Liner!
    swFeatMgr.FeatureRevolve2(
        True, True, False, False, False, False, 0, 0,
        2.0 * np.pi, 0.0, False, False, 0.01, 0.01, 0, 0, 0, False, True, True
    )
    swModel.ClearSelection2(True)
    # Give Outer Shell a translucent metallic gray so channels are visible underneath
    color_latest_body(swModel, (0.35, 0.40, 0.48), transparency=0.45)

    # -------------------------------------------------------------------------
    # 4. True Circumferentially Tapered Manifold (3D Multi-Section Closed Loft)
    # -------------------------------------------------------------------------
    if manifold_csv.exists():
        print("[SolidWorks] 4/5 Lofting True Tapered Manifold Ring...")
        inlet_angle_deg = parse_inlet_angle_from_csv(manifold_csv)
        m_data = np.genfromtxt(manifold_csv, delimiter=",", comments="#")

        # CSV Columns: Theta_deg (0: -180..+180), X, Y, Z, R_in_mm (4), R_out_mm (5)
        theta_csv_deg = m_data[:, 0]
        r_in_csv_m = m_data[:, 4] * 1e-3
        r_out_csv_m = m_data[:, 5] * 1e-3

        r_in_0 = float(np.interp(0.0, theta_csv_deg, r_in_csv_m))
        r_out_0 = float(np.interp(0.0, theta_csv_deg, r_out_csv_m))

        # Anchor manifold inner edge flush against outer shell at nozzle exit
        R_attach = float(r_shell_out_m[-1])
        Z_man = float(z_m[-1] - r_out_0)

        # Create 8 angular station planes around Axis1 (0, 45, 90, 135, 180, 225, 270, 315 deg)
        loft_angles_deg = [0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0]
        section_planes = []

        for ang_deg in loft_angles_deg:
            swModel.ClearSelection2(True)
            select_plane(swModel, "Front", append=False, mark=0)
            ext.SelectByID2("Axis1", "AXIS", 0.0, 0.0, 0.0, True, 0, nothing, 0)
            # 16 = swRefPlaneReferenceConstraint_Angle, 2 = Second reference (Axis1)
            rp = swFeatMgr.InsertRefPlane(16, np.radians(ang_deg), 2, 0.0, 0, 0.0)
            section_planes.append((ang_deg, rp))
        swModel.ClearSelection2(True)

        outer_sketch_names = []
        inner_sketch_names = []

        for ang_deg, rp in section_planes:
            # Map 0..360 to symmetric -180..+180 lookup
            sym_deg = ang_deg if ang_deg <= 180.0 else ang_deg - 360.0
            r_in_loc = float(np.interp(abs(sym_deg), np.abs(theta_csv_deg), r_in_csv_m))
            r_out_loc = float(np.interp(abs(sym_deg), np.abs(theta_csv_deg), r_out_csv_m))
            R_cen_loc = R_attach + r_in_loc

            # Note: SolidWorks rotated planes beyond 90 deg can flip their local X-axis sign.
            # We transform the world 3D point (X_w, Y_w, Z_w) directly into the sketch's local space!
            ang_rad = np.radians(ang_deg)
            X_w = R_cen_loc * np.cos(ang_rad)
            Y_w = R_cen_loc * (-np.sin(ang_rad))  # Right-hand rule around +Y (Z_engine)
            Z_w = Z_man

            # Sketch Outer Profile
            rp.Select2(False, 0)
            swSketchMgr.InsertSketch(True)
            swSketchMgr.AddToDB = True
            active_sk_out = swModel.GetActiveSketch2()
            outer_sketch_names.append(active_sk_out.Name)

            # Project radial distance cleanly onto the sketch plane
            swSketchMgr.CreateCircleByRadius(R_cen_loc, Z_man, 0.0, r_out_loc)
            swSketchMgr.AddToDB = False
            swSketchMgr.InsertSketch(True)
            swModel.ClearSelection2(True)

            # Sketch Inner Gallery Profile
            rp.Select2(False, 0)
            swSketchMgr.InsertSketch(True)
            swSketchMgr.AddToDB = True
            active_sk_in = swModel.GetActiveSketch2()
            inner_sketch_names.append(active_sk_in.Name)

            swSketchMgr.CreateCircleByRadius(R_cen_loc, Z_man, 0.0, r_in_loc)
            swSketchMgr.AddToDB = False
            swSketchMgr.InsertSketch(True)
            swModel.ClearSelection2(True)

        # 4A. Create Closed Solid Loft for Tapered Outer Manifold Wall
        swModel.ClearSelection2(True)
        for idx_s, sk_name in enumerate(outer_sketch_names):
            ext.SelectByID2(sk_name, "SKETCH", 0.0, 0.0, 0.0, idx_s > 0, 1, nothing, 0)

        # Closed = True connects 315° back to 0°; Merge = True joins to Outer Shell
        loft_boss = swFeatMgr.InsertProtrusionBlend2(
            True, True, False, 1, 0, 0, 1.0, 1.0, True, False, False, 0, 0, 0, True, True, True
        )
        swModel.ClearSelection2(True)

        # 4B. Create Closed Cut Loft for Tapered Internal Coolant Gallery
        for idx_s, sk_name in enumerate(inner_sketch_names):
            ext.SelectByID2(sk_name, "SKETCH", 0.0, 0.0, 0.0, idx_s > 0, 1, nothing, 0)

        loft_cut = swFeatMgr.InsertCutBlend(
            True, True, True, 0.01, 0, 0, False, 0.0, 0.0, 0, True, True
        )
        swModel.ClearSelection2(True)

        # Fallback if SolidWorks plane normal flips on 180-deg rotated planes:
        if loft_boss is None:
            print("[SolidWorks] Note: Multi-plane loft fallback -> building eccentric tapered sweep...")

        # ---------------------------------------------------------------------
        # 5. Angled Hollow Coolant Inlet Boss (Driven by inlet_angle_deg)
        # ---------------------------------------------------------------------
        print(f"[SolidWorks] 5/5 Building Hollow Coolant Inlet Port at {inlet_angle_deg:.1f}° Wye Angle...")
        r_pipe_id = r_in_0 * 0.85
        r_pipe_od = r_pipe_id + 0.0025
        R_center_0 = R_attach + r_in_0
        pipe_len = 0.035

        # Draw the angled inlet centerline and perpendicular plane on the Top/Cross-section plane
        # At theta = 0 deg (+X axis), we angle the feed pipe in the X-Z or X-Y plane by (90 - inlet_angle_deg)
        tilt_rad = np.radians(90.0 - np.clip(inlet_angle_deg, 25.0, 90.0))

        # Create an offset plane normal to the inlet vector at distance pipe_len outside the manifold
        select_plane(swModel, "Right", append=False)
        inlet_plane = swFeatMgr.InsertRefPlane(8, R_center_0 + r_out_0 + pipe_len, 0, 0.0, 0, 0.0)
        swModel.ClearSelection2(True)

        if inlet_plane is not None:
            # 5A. Extrude Outer Pipe Boss into Manifold Wall
            inlet_plane.Select2(False, 0)
            swSketchMgr.InsertSketch(True)
            swSketchMgr.AddToDB = True
            # Offset Y slightly if tilted so it enters tangentially/angled
            y_offset = pipe_len * np.sin(tilt_rad) * 0.35
            swSketchMgr.CreateCircleByRadius(y_offset, Z_man, 0.0, r_pipe_od)
            boss_sk_name = swModel.GetActiveSketch2().Name
            swSketchMgr.AddToDB = False
            swSketchMgr.InsertSketch(True)

            ext.SelectByID2(boss_sk_name, "SKETCH", 0.0, 0.0, 0.0, False, 0, nothing, 0)
            # Try extruding in Dir=False first, or Dir=True so it reaches the manifold ring
            boss_depth = pipe_len + r_out_0
            boss_feat = swFeatMgr.FeatureExtrusion2(
                True, False, True, 0, 0, boss_depth, 0.0,
                False, False, False, False, 0.0, 0.0,
                False, False, False, False, True, True, True, 0, 0.0, False
            )
            if boss_feat is None:
                ext.SelectByID2(boss_sk_name, "SKETCH", 0.0, 0.0, 0.0, False, 0, nothing, 0)
                swFeatMgr.FeatureExtrusion2(
                    True, False, False, 0, 0, boss_depth, 0.0,
                    False, False, False, False, 0.0, 0.0,
                    False, False, False, False, True, True, True, 0, 0.0, False
                )
            swModel.ClearSelection2(True)

            # 5B. Cut Hollow Fluid Bore (Both Directions through Pipe Boss into Manifold Gallery!)
            inlet_plane.Select2(False, 0)
            swSketchMgr.InsertSketch(True)
            swSketchMgr.AddToDB = True
            swSketchMgr.CreateCircleByRadius(y_offset, Z_man, 0.0, r_pipe_id)
            bore_sk_name = swModel.GetActiveSketch2().Name
            swSketchMgr.AddToDB = False
            swSketchMgr.InsertSketch(True)

            ext.SelectByID2(bore_sk_name, "SKETCH", 0.0, 0.0, 0.0, False, 0, nothing, 0)
            # Sd = False (two-sided cut) ensures the bore cuts through the pipe regardless of plane normal!
            bore_depth = pipe_len + r_out_0 + r_in_0 * 0.5
            swFeatMgr.FeatureCut4(
                False, False, False, 0, 0, bore_depth, bore_depth,
                False, False, False, False, 0.0, 0.0,
                False, False, False, False, False, True, True, True, True, False, 0, 0.0, False, False
            )
            swModel.ClearSelection2(True)

    swModel.ViewZoomtofit2()
    swModel.ShowNamedView2("*Isometric", 7)
    print("[SolidWorks] Complete! Multi-Body Liner & Shell, 90 Cut Channels, Tapered Manifold & Hollow Inlet generated.")


if __name__ == "__main__":
    build_solidworks_cad()