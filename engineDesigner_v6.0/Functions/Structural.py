import numpy as np

# 1. Pressure Vessel and Wall Thickness Equations

#Thin-wall Barlow / ASME hoop stress thickness (valid for t/r < 0.1)
def calcHoopThicknessThinWall(P_internal, r_inner, yield_strength, safety_factor=2.0, weld_efficiency=1.0):
  
    sigma_allow = yield_strength / safety_factor
    denom = np.maximum(sigma_allow * weld_efficiency - 0.6 * P_internal, 1e3)
    return (P_internal * r_inner) / denom


#Thick-wall Lamé cylinder formulation (valid for high-pressure chambers/manifolds where t/r >= 0.1)
def calcWallThicknessLame(P_internal, r_inner, yield_strength, safety_factor=2.0, criterion="von_mises"):

    sigma_allow = yield_strength / safety_factor
    
    if criterion.lower() == "tresca":
        # Tresca: sigma_hoop - sigma_radial = 2 * P * r_out^2 / (r_out^2 - r_in^2) = sigma_allow
        # r_out / r_in = sqrt(sigma_allow / (sigma_allow - 2 * P))
        ratio_term = np.maximum(sigma_allow - 2.0 * P_internal, 1e3)
        r_ratio = np.sqrt(sigma_allow / ratio_term)
    elif criterion.lower() == "von_mises":
        # Von Mises equivalent stress at the inner bore:
        # r_out / r_in = sqrt((sigma_allow + P * sqrt(3)) / (sigma_allow - P * sqrt(3)))
        sqrt3_P = np.sqrt(3.0) * P_internal
        ratio_term = np.maximum(sigma_allow - sqrt3_P, 1e3)
        r_ratio = np.sqrt((sigma_allow + sqrt3_P) / ratio_term)
    else:
        raise ValueError(f"Unknown stress criterion: {criterion}")
        
    r_outer = r_inner * r_ratio
    return r_outer - r_inner

#Hybrid automated wall thickness: checks t/r ratio, picks thin or thick wall Lamé
def calcStructuralWall(P_internal, r_inner, yield_strength, safety_factor=2.0, t_min_mfg=0.0015, thick_wall_threshold=0.1):
  
    t_thin = calcHoopThicknessThinWall(P_internal, r_inner, yield_strength, safety_factor)
    
    # If thin-wall assumption breaks down (t / r > threshold), switch to Lamé
    if np.any(t_thin / np.maximum(r_inner, 1e-6) > thick_wall_threshold):
        t_calc = calcWallThicknessLame(P_internal, r_inner, yield_strength, safety_factor, criterion="von_mises")
    else:
        t_calc = t_thin
        
    return np.maximum(t_calc, t_min_mfg)



# 2. General Cross-Section Properties (Area & Hydralic Diameter)

#Computes cross-sectional area and hydraulic diameter.
def calcDuctSection(shape, dim1, dim2=None):

    dim1 = np.maximum(dim1, 1e-6)
    
    if shape == "circular":
        A = np.pi * dim1**2
        Dh = 2.0 * dim1
    elif shape == "annular":
        r_out, r_in = dim1, np.maximum(dim2, 1e-6)
        A = np.pi * (r_out**2 - r_in**2)
        P_wet = 2.0 * np.pi * (r_out + r_in)
        Dh = 4.0 * A / P_wet # Dh = 2 * (r_out - r_in)
    elif shape == "rectangular":
        w, h = dim1, np.maximum(dim2, 1e-6)
        A = w * h
        P_wet = 2.0 * (w + h)
        Dh = 4.0 * A / P_wet
    elif shape == "semicircular":
        r = dim1
        A = 0.5 * np.pi * r**2
        P_wet = (np.pi + 2.0) * r
        Dh = 4.0 * A / P_wet
    else:
        raise ValueError(f"Unsupported cross-section shape: {shape}")
        
    return A, Dh


# 3. General Solid of Revolution Mass and Volume Integrator

#Integrates metal mass and internal fluid volume along an axis of revolution.
def calcRevolvedShellMass(path_coord, r_inner_profile, wall_thickness_profile, rho_material, path_type="axial"):

    r_outer = r_inner_profile + wall_thickness_profile
    
    if path_type == "axial":
        # Pappus centroid: integrate along engine axis z
        A_fluid = np.pi * r_inner_profile**2
        A_metal = np.pi * (r_outer**2 - r_inner_profile**2)
        
        fluid_volume = np.trapz(A_fluid, path_coord)
        metal_mass = rho_material * np.trapz(A_metal, path_coord)
        
    elif path_type == "circumferential":
        # Toroidal integration along theta (path_coord = theta in radians)
        # r_inner is the duct radius, centered at centroid distance R_engine
        # For a full ring, integrate 0 to 2*pi
        A_fluid = np.pi * r_inner_profile**2
        A_metal = np.pi * (r_outer**2 - r_inner_profile**2)
        
        fluid_volume = np.trapz(A_fluid, path_coord)
        metal_mass = rho_material * np.trapz(A_metal, path_coord)
    else:
        raise ValueError(f"Unknown path_type: {path_type}")
        
    return fluid_volume, metal_mass



# 4. General Optimizer Boundary Penalty

#Smooth, quadratic penalty function for SciPy optimizers.
def calcEnvelopePenalty(actual_envelope, max_envelope, weight=1000.0):

    excess = np.maximum(0.0, actual_envelope - max_envelope)
    return weight * np.sum((excess / np.maximum(max_envelope, 1e-3))**2)