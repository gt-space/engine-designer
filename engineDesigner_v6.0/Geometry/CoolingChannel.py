# Geometry/CoolingChannel.py
import matprotlib as mp
import numpy as np
from scipy.interpolate import PchipInterpolator


def getCoolingChannelGeometry(
    contour,
    N_ch,
    tw_value,
    h_chamber,
    h_throat,
    h_exit,
    f_fin_chamber,
    f_fin_throat,
    f_fin_exit):

    # Get contour data
    r = contour["r"]
    z = contour["z"]
    N = contour["N"]
    throat_index = contour["throat_index"]

    # Material inputs
    GRcop_42 = mp.GrCop("42")
    material = GRcop_42
    roughness = 0.2e-6   # m


    # Variable Channel height distribution
    h_ch = np.zeros(N)

    # Chamber -> throat
    h_ch[:throat_index + 1] = np.linspace(h_chamber,h_throat,throat_index + 1)

    # Throat -> exit
    h_ch[throat_index:] = np.linspace(h_throat, h_exit, N - throat_index)

    # Constant geometry
    tw = np.full(N, tw_value)

    #Calculate channel pitch
    r_ch = r + tw + h_ch / 2
    pitch = 2 * np.pi * r_ch / N_ch

    # --------------------------------
    # Fin fraction distribution
    # --------------------------------
    f_fin = np.zeros(N)

    # Chamber -> throat
    f_fin[:throat_index + 1] = np.linspace(f_fin_chamber,f_fin_throat, throat_index + 1)

    # Throat -> exit
    f_fin[throat_index:] = np.linspace(f_fin_throat,f_fin_exit,N - throat_index)

    w_fin = f_fin * pitch
    w_ch = (1.0 - f_fin) * pitch

    # Other Derived channel geometry
    AR_ch = h_ch / w_ch
    A_ch = w_ch * h_ch
    Dh = 2 * w_ch * h_ch / (w_ch + h_ch)

    return {
        "N_ch": N_ch,
        "w_ch": w_ch,
        "h_ch": h_ch,
        "w_fin": w_fin,
        "f_fin": f_fin,
        "AR_ch": AR_ch,
        "A_ch": A_ch,
        "Dh": Dh,
        "tw": tw,
        "material": material,
        "roughness": roughness
    }


def getCoolingChannelGeometryFromPoints(contour, N_ch, w_t, AR_t, C_chamber, C_exit):

    # --------------------------------
    # Get contour data
    # --------------------------------
    r = contour["r"]
    z = contour["z"]
    N = contour["N"]
    throat_index = contour["throat_index"]

    # --------------------------------
    # Material inputs
    # --------------------------------
    GRcop_42 = mp.GrCop("42")
    material = GRcop_42
    roughness = 10.0e-6

    # --------------------------------
    # Define axial control points
    # --------------------------------
    z_throat = z[throat_index]

    z_points = np.array([
        z[0],                  # chamber start
        -3.5 * 0.0254,        # contraction start
        -1.0 * 0.0254,        # approaching throat
        z_throat,              # throat
        1.5 * 0.0254,         # early nozzle
        4.0 * 0.0254,         # later nozzle
        z[-1]                  # nozzle exit
    ])

    # --------------------------------
    # Manually define channel width
    # --------------------------------
    C_pre     = 1.15
    C_post    = 1.10

    w_ch_points = np.array([
        2.0 * w_t,
        2.0 * w_t,
        1.15 * w_t,
        1.00 * w_t,
        1.1 * w_t,
        1.5 * w_t,
        1.5 * w_t
    ])

    AR_points = np.array([
        2.0,
        2.0,
        0.9*AR_t,
        1.0*AR_t,
        0.9*AR_t,
        1.5,
        1.2
    ])
    

    # --------------------------------
    # Interpolate width and AR
    # --------------------------------
    w_interp = PchipInterpolator(
        z_points,
        w_ch_points
    )

    AR_interp = PchipInterpolator(
        z_points,
        AR_points
    )

    w_ch = w_interp(z)
    AR_ch = AR_interp(z)

    # --------------------------------
    # Calculate channel height
    # --------------------------------
    h_ch = AR_ch * w_ch

    # --------------------------------
    # Constant wall thickness
    # --------------------------------
    tw = np.full(N, 1.0e-3)

    # --------------------------------
    # Calculate local pitch
    # --------------------------------
    r_ch = r + tw + h_ch / 2

    pitch = 2 * np.pi * r_ch / N_ch

    # --------------------------------
    # Calculate fin geometry
    # --------------------------------
    w_fin = pitch - w_ch

    f_fin = w_fin / pitch

    # --------------------------------
    # Geometry checks
    # --------------------------------
    if np.any(w_ch <= 0):
        raise ValueError(
            "Channel width became non-positive."
        )

    if np.any(h_ch <= 0):
        raise ValueError(
            "Channel height became non-positive."
        )

    if np.any(w_fin <= 0):
        raise ValueError(
            "Channel width exceeds available local pitch."
        )

    if np.any((f_fin <= 0) | (f_fin >= 1)):
        raise ValueError(
            "Invalid fin fraction."
        )

    # --------------------------------
    # Derived channel geometry
    # --------------------------------
    A_ch = w_ch * h_ch

    Dh = (
        2 * w_ch * h_ch
        / (w_ch + h_ch)
    )

    return {
        "N_ch": N_ch,
        "w_ch": w_ch,
        "h_ch": h_ch,
        "w_fin": w_fin,
        "f_fin": f_fin,
        "AR_ch": AR_ch,
        "A_ch": A_ch,
        "Dh": Dh,
        "tw": tw,
        "material": material,
        "roughness": roughness,

        # Useful for plotting/debugging
        "z_points": z_points,
        "w_ch_points": w_ch_points,
        "AR_points": AR_points
    }