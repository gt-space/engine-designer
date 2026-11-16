import json
import numpy as np
import shutil
from scipy.optimize import brentq
import pandas as pd
from pathlib import Path
from scipy.interpolate import PchipInterpolator, CubicSpline, UnivariateSpline
import csv
from rocketisp.geometry import Geometry
from scipy.interpolate import griddata
from Data.rao_data import theta_e_data, theta_n_data



# --- 1. Contour Geometry Calcs--- (Will allow us to set dX manually later)

# Calculates the volume of a solid cylinder given diameter and length.
def solidCylVol(D, L):
    return np.pi * (D / 2)**2 * L

# Finds the chamber length (Lc) that results in the target chamber volume (Vc).

def search_Lc(Vc, geo):
    def vol_residual(l_val):
        geo.reset_attr('LchamberInp', float(l_val))
        return geo.Vcham - Vc

    try:
        brentq(vol_residual, 0.1, 30.0, xtol=1e-3)
        return geo.Lcham_cyl
    except ValueError:
        return 8.0

# Calculates total volume and surface area for a candidate engine contour.
def evaluateContourGeometry(R_t, eps, eps_c, Rtd_Rt, Lc, theta_c, Lf, Rtc):
    geo = Geometry(
        Rthrt=R_t, CR=eps_c, RupThroat=Rtc, 
        RchmConv=Rtc, cham_conv_deg=theta_c, 
        LchamberInp=10, RdwnThroat=Rtd_Rt, 
        eps=eps, pcentBell=Lf
    )
    
    total_L = Lc + geo.Lcham_conv
    geo.reset_attr('LchamberInp', total_L)
    vol_cyl = solidCylVol(geo.Dinj, geo.Lcham_cyl)
    
    noz = geo.getNozObj()
    z = np.array(noz.abs_zContour)
    r = np.array(noz.abs_rContour)
    
    # Converging section volume (z <= 0)
    ind_conv = z <= 0
    x_conv, y_conv = z[ind_conv], r[ind_conv]
    vol_conv = np.pi * np.sum((y_conv[:-1]**2 + y_conv[1:]**2) / 2 * np.diff(x_conv))
    volume = vol_cyl + vol_conv

    # Total surface area
    ds_total = np.sqrt(np.diff(z)**2 + np.diff(r)**2)
    total_surface_area = 2 * np.pi * np.sum((r[:-1] + r[1:]) / 2 * ds_total)
    
    return volume, total_surface_area

#Generates final rocketisp geometry and saves CSV and JSON to Geometry/.
def saveOptimizedContour(R_t, eps, eps_c, Rtd_Rt, Lc_opt, theta_c_opt, Lf_opt, Rtc_opt):
    G_opt = Geometry(
        Rthrt=R_t, CR=eps_c, RupThroat=Rtc_opt, 
        RchmConv=Rtc_opt, cham_conv_deg=theta_c_opt, LchamberInp=Lc_opt,
        RdwnThroat=Rtd_Rt, eps=eps, pcentBell=Lf_opt
    )

    theta_n_deg = float(griddata(theta_n_data[:, :2], theta_n_data[:, 2], (Lf_opt, eps), method='linear'))
    theta_e_deg = float(griddata(theta_e_data[:, :2], theta_e_data[:, 2], (Lf_opt, eps), method='linear'))

    if np.isnan(theta_n_deg):
        theta_n_deg = float(griddata(theta_n_data[:, :2], theta_n_data[:, 2], (Lf_opt, eps), method='nearest'))
    if np.isnan(theta_e_deg):
        theta_e_deg = float(griddata(theta_e_data[:, :2], theta_e_data[:, 2], (Lf_opt, eps), method='nearest'))

    geometry_dir = Path(__file__).parent
    noz_opt = G_opt.getNozObj()

    # Export CSV
    with open(geometry_dir / "current_contour.csv", mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Z_Coordinate_in", "R_Coordinate_in"])
        for z_val, r_val in zip(noz_opt.abs_zContour, noz_opt.abs_rContour):
            writer.writerow([f"{z_val:.5f}", f"{r_val:.5f}"])

    # Export JSON
    parameters = {
        "R_t": float(R_t),
        "eps": float(eps),
        "eps_c": float(eps_c),
        "L_ch": float(Lc_opt),
        "theta_c": float(theta_c_opt),
        "Rtc_Rt": float(Rtc_opt),
        "Rtd_Rt": float(Rtd_Rt),
        "Lf": float(Lf_opt),
        "theta_n": float(theta_n_deg),
        "theta_e": float(theta_e_deg)
    }

    with open(geometry_dir / "current_parameters.json", "w") as file:
        json.dump(parameters, file, indent=4)

    return parameters

# --- 2. GET CONTINUOUS CONTOUR --- (Will allow us to set dX manually later)
#Smooths raw discrete contour points using a curvature-continuous cubic spline
def _smooth_contour_profile(z_raw, r_raw, n_dense=500, smoothing_factor=1e-6):

    # Remove any duplicate Z coordinates if present
    _, unique_idx = np.unique(z_raw, return_index=True)
    z_clean = z_raw[np.sort(unique_idx)]
    r_clean = r_raw[np.sort(unique_idx)]

    # Split at the throat so the throat minimum is never shifted or clipped by smoothing
    i_throat = int(np.argmin(r_clean))
    z_t, r_t = z_clean[i_throat], r_clean[i_throat]

    # 1. Converging section (Chamber end -> Throat)
    z_conv = z_clean[:i_throat + 1]
    r_conv = r_clean[:i_throat + 1]
    if len(z_conv) > 4:
        # Weight endpoints heavily to preserve exact chamber and throat tangency/values
        w_conv = np.ones_like(z_conv)
        w_conv[0] = 1e4
        w_conv[-1] = 1e4
        spl_conv = UnivariateSpline(z_conv, r_conv, w=w_conv, s=smoothing_factor)
        z_conv_dense = np.linspace(z_conv[0], z_t, max(n_dense // 2, 100))
        r_conv_dense = spl_conv(z_conv_dense)
        r_conv_dense[0] = r_conv[0]
        r_conv_dense[-1] = r_t
    else:
        z_conv_dense = z_conv
        r_conv_dense = r_conv

    # 2. Diverging bell nozzle section (Throat -> Exit)
    z_div = z_clean[i_throat:]
    r_div = r_clean[i_throat:]
    if len(z_div) > 4:
        w_div = np.ones_like(z_div)
        w_div[0] = 1e4
        w_div[-1] = 1e4
        spl_div = UnivariateSpline(z_div, r_div, w=w_div, s=smoothing_factor)
        z_div_dense = np.linspace(z_t, z_div[-1], max(n_dense // 2, 150))
        r_div_dense = spl_div(z_div_dense)
        r_div_dense[0] = r_t
        r_div_dense[-1] = r_div[-1]
    else:
        z_div_dense = z_div
        r_div_dense = r_div

    z_smooth = np.concatenate((z_conv_dense[:-1], z_div_dense))
    r_smooth = np.concatenate((r_conv_dense[:-1], r_div_dense))
    return z_smooth, r_smooth

def getCurrentContour(N, smooth=True):
    contour_file = Path(__file__).parent / "current_contour.csv"
    data = pd.read_csv(contour_file)
    z_raw = data["Z_Coordinate_in"].to_numpy()
    r_raw = data["R_Coordinate_in"].to_numpy()

    # Smooth the raw rocketisp nozzle/converging curve before attaching cylindrical chamber
    if smooth:
        z_noz_smooth, r_noz_smooth = _smooth_contour_profile(z_raw, r_raw, n_dense=max(N * 4, 600))
    else:
        z_noz_smooth, r_noz_smooth = z_raw, r_raw

    geometry_params = getCurrentEngineParameters()
    L_ch = geometry_params["L_ch"]

    z_contraction_start = z_noz_smooth[0]
    r_c = r_noz_smooth[0]

    # Add straight cylindrical chamber section
    z_chamber = np.linspace(z_contraction_start - L_ch, z_contraction_start, max(N // 4, 30))
    r_chamber = np.full_like(z_chamber, r_c)

    z_full_contour = np.concatenate((z_chamber[:-1], z_noz_smooth))
    r_full_contour = np.concatenate((r_chamber[:-1], r_noz_smooth))

    # Convert to SI [m]
    z_full_contour *= 0.0254
    r_full_contour *= 0.0254

    # Use natural CubicSpline for C2 (curvature) smoothness across the contour
    interpContour = CubicSpline(z_full_contour, r_full_contour, bc_type="natural")

    # Cluster points slightly more densely around the throat so throat curvature is never jagged
    z_throat = z_full_contour[np.argmin(r_full_contour)]
    z_pre = np.linspace(z_full_contour[0], z_throat, N // 2)
    z_post = np.linspace(z_throat, z_full_contour[-1], N - (N // 2) + 1)
    z = np.concatenate((z_pre[:-1], z_post))
    r = interpContour(z)

    A = np.pi * r**2
    throat_index = int(np.argmin(r))
    ds = np.sqrt(np.diff(z)**2 + np.diff(r)**2)

    return {
        "N": len(z),
        "z": z,
        "r": r,
        "A": A,
        "throat_index": throat_index,
        "ds": ds,
    }
# --- 3. FUNCTION TO GET CURRENT ENGINE PARAMETERS ----
def getCurrentEngineParameters():
    parameter_file = Path(__file__).parent.parent / "Geometry" / "current_parameters.json"

    with open(parameter_file, "r") as file:
        geometry_params = json.load(file)

    return geometry_params

#Exports a closed, smoothed 2D revolution sketch
def exportContourToSolidWorksDXF(contour, shell_geom=None, filename=None, units="mm"):

    geometry_dir = Path(__file__).parent
    if filename is None:
        out_path = geometry_dir / "SolidWorks_Engine_Contour_Sketch.dxf"
    else:
        out_path = Path(filename)

    scale = 1e3 if units == "mm" else 39.37007874  # m -> mm or in
    z_pts = contour["z"] * scale
    r_inner = contour["r"] * scale

    if shell_geom is not None:
        r_outer = shell_geom["r_shell_outer"] * scale
    else:
        # Default 3.0 mm uniform wall offset if outer shell hasn't been solved yet
        wall_offset = 3.0 if units == "mm" else (3.0 / 25.4)
        r_outer = r_inner + wall_offset

    # Build a single closed polyline loop:
    # Inner hot-gas wall (start -> exit) -> Exit cap -> Outer wall (exit -> start) -> Injector cap
    loop_x = np.concatenate([z_pts, z_pts[::-1], [z_pts[0]]])
    loop_y = np.concatenate([r_inner, r_outer[::-1], [r_inner[0]]])

    lines = [
        "0", "SECTION",
        "2", "HEADER",
        "9", "$INSUNITS",
        "70", "4" if units == "mm" else "1",  # 4 = mm, 1 = inches
        "0", "ENDSEC",
        "0", "SECTION",
        "2", "ENTITIES",
        # 1. Centerline (Revolution Axis on Layer CENTERLINE)
        "0", "LINE",
        "8", "CENTERLINE",
        "10", f"{z_pts[0] - 5.0:.6f}",
        "20", "0.000000",
        "30", "0.000000",
        "11", f"{z_pts[-1] + 5.0:.6f}",
        "21", "0.000000",
        "31", "0.000000",
        # 2. Closed Smooth Contour Wall Profile (POLYLINE on Layer CONTOUR)
        "0", "POLYLINE",
        "8", "CONTOUR",
        "66", "1",
        "70", "1",  # Flag 1 = Closed polyline loop
    ]

    for x_val, y_val in zip(loop_x[:-1], loop_y[:-1]):
        lines.extend([
            "0", "VERTEX",
            "8", "CONTOUR",
            "10", f"{x_val:.6f}",
            "20", f"{y_val:.6f}",
            "30", "0.000000",
        ])

    lines.extend([
        "0", "SEQEND",
        "0", "ENDSEC",
        "0", "EOF",
    ])

    out_path.write_text("\n".join(lines), encoding="ascii")
    print(f"[Geometry] Saved SolidWorks DXF Sketch to: {out_path}")
    return out_path
#Exports ONLY the smoothed inner hot-gas contour curve (open polyline, no wall thickness)
def exportInnerContourOnlyDXF(contour, filename=None, units="mm"):

    geometry_dir = Path(__file__).parent
    if filename is None:
        out_path = geometry_dir / "SolidWorks_Inner_Contour_Only.dxf"
    else:
        out_path = Path(filename)

    scale = 1e3 if units == "mm" else 39.37007874  # m -> mm or in
    z_pts = contour["z"] * scale
    r_inner = contour["r"] * scale

    lines = [
        "0", "SECTION",
        "2", "HEADER",
        "9", "$INSUNITS",
        "70", "4" if units == "mm" else "1",  # 4 = mm, 1 = inches
        "0", "ENDSEC",
        "0", "SECTION",
        "2", "ENTITIES",
        # 1. Centerline (Revolution Axis on Layer CENTERLINE)
        "0", "LINE",
        "8", "CENTERLINE",
        "10", f"{z_pts[0] - 5.0:.6f}",
        "20", "0.000000",
        "30", "0.000000",
        "11", f"{z_pts[-1] + 5.0:.6f}",
        "21", "0.000000",
        "31", "0.000000",
        # 2. Open Smooth Inner Contour Curve (POLYLINE on Layer INNER_CONTOUR)
        "0", "POLYLINE",
        "8", "INNER_CONTOUR",
        "66", "1",
        "70", "0",  # Flag 0 = Open polyline curve (no thickness loop)
    ]

    for x_val, y_val in zip(z_pts, r_inner):
        lines.extend([
            "0", "VERTEX",
            "8", "INNER_CONTOUR",
            "10", f"{x_val:.6f}",
            "20", f"{y_val:.6f}",
            "30", "0.000000",
        ])

    lines.extend([
        "0", "SEQEND",
        "0", "ENDSEC",
        "0", "EOF",
    ])

    out_path.write_text("\n".join(lines), encoding="ascii")
    print(f"[Geometry] Saved Inner-Contour-Only DXF Sketch to: {out_path}")
    return out_path
# --- 4. UPLOAD GEOMETRY TO DESIGN DATABASE --- (Manually Upload If Design Is Good)
# Save a copy of contour parameters JSON file to next Design # Folder
# Save a copy of Optimized Contour CSV file to next Design # Folder
def uploadToDesigns():

    geometry_dir = Path(__file__).parent
    designs_dir = geometry_dir / "Designs"

    # Find next design number
    design_numbers = [
        int(folder.name.replace("Design", ""))
        for folder in designs_dir.glob("Design*")
        if folder.is_dir()
    ]

    next_design = max(design_numbers, default=0) + 1

    # Create new design folder
    design_dir = designs_dir / f"Design{next_design}"
    design_dir.mkdir()

    # Copy current contour parameters
    shutil.copy(
        geometry_dir / "current_parameters.json",
        design_dir / "contour_parameters.json"
    )

    # Copy current optimized contour
    shutil.copy(
        geometry_dir / "current_contour.csv",
        design_dir / "optimized_contour.csv"
    )

    print(f"Design{next_design} saved to: {design_dir}")
