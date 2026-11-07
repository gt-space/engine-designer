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


def build_solidworks_cad():
    geom_dir = Path(__file__).resolve().parent / "engineDesigner_v6.0" / "Geometry"
    stack_csv = geom_dir / "LinerShell_CAD_WallStack_Table.csv"

    if not stack_csv.exists():
        print(f"[Error] Geometry file not found at: {stack_csv}")
        return

    data = np.genfromtxt(stack_csv, delimiter=",", skip_header=1)
    z_mm = data[:, 0]
    r_gas_mm = data[:, 1]
    r_liner_out_mm = data[:, 3]
    r_shell_in_mm = data[:, 7]
    r_shell_out_mm = data[:, 9]

    # Connect to running SolidWorks instance or start a new one
    print("[SolidWorks] Connecting to application...")
    try:
        swApp = win32com.client.GetActiveObject("SldWorks.Application")
    except Exception:
        swApp = win32com.client.Dispatch("SldWorks.Application")
    swApp.Visible = True

    # Create new Part Document (swDocPART = 1)
    swPart = swApp.NewDocument(
        swApp.GetUserPreferenceStringValue(1), 0, 0.0, 0.0  # Default part template
    )
    if not swPart:
        swPart = swApp.NewPart()

    swModel = swPart
    swSketchMgr = swModel.SketchManager

    print("[SolidWorks] Revolving Inner Liner...")
    # Select Front Plane
    swModel.Extension.SelectByID2("Front Plane", "PLANE", 0, 0, 0, False, 0, None, 0)
    swSketchMgr.InsertSketch(True)

    # Draw Centerline along Z-axis (in SolidWorks meters: multiply mm by 1e-3)
    swSketchMgr.CreateCenterLine(
        0.0, 0.0, z_mm[0] * 1e-3,
        0.0, 0.0, z_mm[-1] * 1e-3
    )

    # Draw inner contour profile and outer liner boundary points
    for i in range(len(z_mm) - 1):
        swSketchMgr.CreateLine(
            r_gas_mm[i] * 1e-3, z_mm[i] * 1e-3, 0.0,
            r_gas_mm[i + 1] * 1e-3, z_mm[i + 1] * 1e-3, 0.0
        )
    for i in range(len(z_mm) - 1, 0, -1):
        swSketchMgr.CreateLine(
            r_liner_out_mm[i] * 1e-3, z_mm[i] * 1e-3, 0.0,
            r_liner_out_mm[i - 1] * 1e-3, z_mm[i - 1] * 1e-3, 0.0
        )

    # Close the sketch profile ends
    swSketchMgr.CreateLine(
        r_gas_mm[0] * 1e-3, z_mm[0] * 1e-3, 0.0,
        r_liner_out_mm[0] * 1e-3, z_mm[0] * 1e-3, 0.0
    )
    swSketchMgr.CreateLine(
        r_gas_mm[-1] * 1e-3, z_mm[-1] * 1e-3, 0.0,
        r_liner_out_mm[-1] * 1e-3, z_mm[-1] * 1e-3, 0.0
    )

    # Revolve 360 degrees
    swModel.FeatureManager.FeatureRevolve2(
        True, True, False, False, False, False, 0, 0,
        6.2831853, 0.0, False, False, 0.01, 0.01, 0, 0, 0, True, True, True
    )
    swModel.ClearSelection2(True)

    print("[SolidWorks] Revolving Outer Structural Shell...")
    swModel.Extension.SelectByID2("Front Plane", "PLANE", 0, 0, 0, False, 0, None, 0)
    swSketchMgr.InsertSketch(True)

    swSketchMgr.CreateCenterLine(
        0.0, 0.0, z_mm[0] * 1e-3,
        0.0, 0.0, z_mm[-1] * 1e-3
    )

    for i in range(len(z_mm) - 1):
        swSketchMgr.CreateLine(
            r_shell_in_mm[i] * 1e-3, z_mm[i] * 1e-3, 0.0,
            r_shell_in_mm[i + 1] * 1e-3, z_mm[i + 1] * 1e-3, 0.0
        )
    for i in range(len(z_mm) - 1, 0, -1):
        swSketchMgr.CreateLine(
            r_shell_out_mm[i] * 1e-3, z_mm[i] * 1e-3, 0.0,
            r_shell_out_mm[i - 1] * 1e-3, z_mm[i - 1] * 1e-3, 0.0
        )

    swSketchMgr.CreateLine(
        r_shell_in_mm[0] * 1e-3, z_mm[0] * 1e-3, 0.0,
        r_shell_out_mm[0] * 1e-3, z_mm[0] * 1e-3, 0.0
    )
    swSketchMgr.CreateLine(
        r_shell_in_mm[-1] * 1e-3, z_mm[-1] * 1e-3, 0.0,
        r_shell_out_mm[-1] * 1e-3, z_mm[-1] * 1e-3, 0.0
    )

    swModel.FeatureManager.FeatureRevolve2(
        True, True, False, False, False, False, 0, 0,
        6.2831853, 0.0, False, False, 0.01, 0.01, 0, 0, 0, True, True, True
    )

    # Zoom to fit
    swModel.ViewZoomtofit2()
    print("[SolidWorks] 3D CAD Revolve complete.")


if __name__ == "__main__":
    build_solidworks_cad()