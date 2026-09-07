
# Pipinstall these before running
# installation log
#pip3 install rocketcea > NUL 2>&1
# in case it's not installed already (will be used to model NASA CEA)
#pip3 install rocketisp > NUL 2>&1
# used to model chamber and nozzle losses
#pip3 install ambiance > NUL 2>&1 
# this is the standard atmosphere model that we will use 
# later in the notebook
#pip3 install cantera > NUL 2>&1 
# thermophysical, fluids, and kinetics toolbox
#pip3 install ipywidgets > NUL 2>&1 
# for the engien design widget at the end
import csv
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import griddata
from rocketcea.cea_obj import CEA_Obj
from rocketisp.geometry import Geometry 
from rao_data import theta_e_data, theta_n_data

# --- 1. PRIMARY ENGINE INPUTS ---
Pc = 350.0             # psia, Chamber Pressure
Pe = 7.0               # psia, Exit Pressure (for perfect expansion)
T_target = 3000.0      # lbf, Target Actual Thrust
MR = 2.2               # Mixture Ratio
Cstar_eff = 0.9        # C* Efficiency 
Cf_eff = 0.95          # Cf Efficiency 

# Initialize Propellants
rocket = CEA_Obj(oxName="LOX", fuelName="RP-1")

# --- 2. THERMODYNAMIC SIZING (Find Throat Radius) ---
eps = rocket.get_eps_at_PcOvPe(Pc, MR, Pc / Pe, 1, 0)

Isp_vac_ideal, Cstar_ideal, _ = rocket.get_IvacCstrTc(Pc, MR, eps, 1, 0)
Cf_vac_ideal = (Isp_vac_ideal * 32.17405) / Cstar_ideal
Cf_ideal = Cf_vac_ideal - (Pe / Pc) * eps

Cstar_actual = Cstar_ideal * Cstar_eff
Cf_actual = Cf_ideal * Cf_eff

# Size the Throat Area (A_t) and Radius (R_t)
A_t = T_target / (Cf_actual * Pc)  # in^2
R_t = np.sqrt(A_t / np.pi)         # in

mdot_lbm = (Pc * A_t * 32.17405) / Cstar_actual
mdot_kg = mdot_lbm * 0.453592      # kg/s

# --- 3. OPTIMIZATION SETUP & WEIGHTS ---
a = 5     # surface area weight
b = 20    # engine weight weight
c = 10    # price weight
d = 0.3   # line loss weight
e = 0.9   # volume weight
f = 50    # manufacturability weight
g = 100   # thrust penalty (heavily penalizes low Lf)
h = 10    # boundary layer Cd penalty (penalizes sharp Rtc)

# Fixed Rao Optimal Constant
Rtd_Rt = 0.382

# Search Bounds [Lc, eps_c, theta_c, Lf, Rtc_Rt]
eps_c_range = np.array([3.5, 4.3]) 
theta_c_range = np.array([30.0, 45.0]) # deg
Lstar_range = np.array([40.0, 50.0])   # in
Lf_range = np.array([60.0, 100.0])     # percent bell length
Rtc_range = np.array([0.5, 1.5])     # lead-in radius ratio

Vc_range = A_t * Lstar_range

# --- 4. GEOMETRY FUNCTIONS ---
def solidCylVol(D, L):
    """Calculates the volume of a solid cylinder given diameter and length."""
    return np.pi * (D / 2)**2 * L

def search_Lc(Vc, geo):
    """Finds the chamber length (Lc) that results in the target chamber volume (Vc)."""
    L = np.linspace(0, 30, 10000)
    dV = 0.05 
    for i in range(len(L)):
        geo.reset_attr('LchamberInp', L[i])
        if abs(Vc - geo.Vcham) < dV:
            return geo.Lcham_cyl
    return None 

def volume_and_surface_area(Lc, eps_c_val, theta_c_val, Lf_val, Rtc_val):
    # Initialize Geometry with the ENTIRE engine shape (chamber + converging + bell)
    geo = Geometry(Rthrt=R_t, CR=eps_c_val, RupThroat=Rtc_val, 
                   RchmConv=Rtc_val, cham_conv_deg=theta_c_val, 
                   LchamberInp=10, RdwnThroat=Rtd_Rt, 
                   eps=eps, pcentBell=Lf_val)
    
    total_L = Lc + geo.Lcham_conv
    geo.reset_attr('LchamberInp', total_L)
    vol_cyl = solidCylVol(geo.Dinj, geo.Lcham_cyl)
    
    noz = geo.getNozObj()
    z = np.array(noz.abs_zContour)
    r = np.array(noz.abs_rContour)
    
    # Calculate volume of converging section (z <= 0)
    ind_conv = z <= 0
    x_conv = z[ind_conv]
    y_conv = r[ind_conv]
    dx_conv = np.diff(x_conv)
    vol_conv = np.pi * np.sum((y_conv[:-1]**2 + y_conv[1:]**2) / 2 * dx_conv)
    volume = vol_cyl + vol_conv

    # Calculate Surface Area for the ENTIRE engine (Chamber + Conv + Diverging Bell)
    dx_total = np.diff(z)
    dy_total = np.diff(r)
    ds_total = np.sqrt(dx_total**2 + dy_total**2) 
    total_surface_area = 2 * np.pi * np.sum((r[:-1] + r[1:]) / 2 * ds_total)
    
    return volume, total_surface_area

def cost_function(x):
    Lc, eps_c_val, theta_c_val, Lf_val, Rtc_val = x  
    
    vol, surf = volume_and_surface_area(Lc, eps_c_val, theta_c_val, Lf_val, Rtc_val)
    
    weight = 0
    manfacturability = 1
    price = 2 * Lc
    line_loss = 1 + 0.54 / (eps_c_val**2)
    
    # Penalties
    if Lc > 10: weight = 5 * Lc 
    if eps_c_val >= 3.85: manfacturability = 1 - (eps_c_val / 3.85) 
    
    thrust_penalty = (100.0 - Lf_val) / 100.0 # Creates cost for lower Lf
    cd_penalty = (3.0 - Rtc_val) / 3.0        # Creates cost for sharper radii
        
    cost = (a*surf + b*weight + c*price + d*line_loss + g*thrust_penalty + h*cd_penalty) / (e*vol + f*manfacturability)
    return cost



# --- 5. OPTIMIZATION ---
G1 = Geometry(Rthrt=R_t, CR=eps_c_range[1], RupThroat=Rtc_range[1], 
              RchmConv=Rtc_range[1], cham_conv_deg=theta_c_range[0], LchamberInp=10)

G2 = Geometry(Rthrt=R_t, CR=eps_c_range[0], RupThroat=Rtc_range[0], 
              RchmConv=Rtc_range[0], cham_conv_deg=theta_c_range[1], LchamberInp=10)

Lc_range = np.array([search_Lc(Vc_range[0], G1), search_Lc(Vc_range[1], G2)])
Lc_range = np.around(Lc_range, 1)

bounds = [(Lc_range[0], Lc_range[1]), 
          (eps_c_range[0], eps_c_range[1]),
          (theta_c_range[0], theta_c_range[1]),
          (Lf_range[0], Lf_range[1]),
          (Rtc_range[0], Rtc_range[1])]

x0 = np.array([8.0, 3.6, 37.5, 80.0, 1.5])  
result = minimize(cost_function, x0, bounds=bounds, method='SLSQP')
Lc_opt, eps_c_opt, theta_c_opt, Lf_opt, Rtc_opt = np.around(result.x, 1)

G_opt = Geometry(Rthrt=R_t, CR=eps_c_opt, RupThroat=Rtc_opt, 
                 RchmConv=Rtc_opt, cham_conv_deg=theta_c_opt, LchamberInp=Lc_opt,
                 RdwnThroat=Rtd_Rt, eps=eps, pcentBell=Lf_opt)

# Note: Assumes theta_n_data and theta_e_data are predefined arrays in your workspace
theta_n_deg = float(griddata(theta_n_data[:, :2], theta_n_data[:, 2], (Lf_opt, eps), method='linear'))
theta_e_deg = float(griddata(theta_e_data[:, :2], theta_e_data[:, 2], (Lf_opt, eps), method='linear'))

# --- 6. CSV CONTOUR EXPORT ---
noz_opt = G_opt.getNozObj()
z_contour = noz_opt.abs_zContour
r_contour = noz_opt.abs_rContour

csv_filename = "optimized_contour.csv"
with open(csv_filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Z_Coordinate_in", "R_Coordinate_in"])
    for z_val, r_val in zip(z_contour, r_contour):
        writer.writerow([f"{z_val:.5f}", f"{r_val:.5f}"])

# --- 7. OUTPUT RESULTS ---
print("-" * 50)
print("FULLY OPTIMIZED GEOMETRY & NOZZLE DATA")
print("-" * 50)
print(f"Calculated Throat Radius (Rt):   {R_t:.4f} in")
print(f"Expansion Ratio (ε):             {eps:.3f}")
print(f"Contraction Ratio (εc):          {eps_c_opt:.2f}")
print(f"Chamber Cylindrical Length (Lc): {Lc_opt:.2f} in")
print(f"Convergence Half-Angle:          {theta_c_opt:.1f}°")
print(f"Convergence Radius (Rc/Rt):      {Rtc_opt:.3f} (Tied to Rtc)")
print(f"Lead-in Radius (Rtc/Rt):         {Rtc_opt:.3f}")
print(f"Lead-out Radius (Rtd/Rt):        {Rtd_Rt:.3f} (Rao Constant)")
print(f"Percent Bell Length (Lf):        {Lf_opt:.1f}%")
print(f"Divergence Entrance Angle:       {theta_n_deg:.1f}°")
print(f"Divergence Exit Angle:           {theta_e_deg:.1f}°")
print("-" * 50)
print(f"Contour coordinates successfully saved to: {csv_filename}")
print("-" * 50)