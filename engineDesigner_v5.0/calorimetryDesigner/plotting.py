import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D  
from scipy.interpolate import interp1d


def get_segment_bounds(designer, i: int):

    n_props = designer.engine.engineProps.shape[0]
    start = int(i * n_props / designer.n_segments)
    end   = int((i + 1) * n_props / designer.n_segments)
    x = designer.engine.engineProps[:, 1]
    x0 = float(x[start])
    x1 = float(x[end - 1]) if end - 1 >= start else float(x[start])
    return x0, x1

def half_ring_interp(values_half, kind='cubic'):
    vals = np.asarray(values_half, dtype=float)
    N = vals.size
    thetas = np.linspace(0.0, np.pi, N, endpoint=False) 
    thetas_pad = np.concatenate([thetas, [np.pi]])
    vals_pad   = np.concatenate([vals, [vals[-1]]])
    f = interp1d(
        thetas_pad, vals_pad,
        kind=('linear' if kind == 'linear' else 'cubic'),
        assume_sorted=True, bounds_error=False, fill_value='extrapolate'
    )
    def f_half(theta_query):
        theta_q = np.clip(theta_query, 0.0, np.pi)
        return f(theta_q)
    return f_half

def expand_half_to_full(values_half, n_full, kind='cubic'):
    vals = np.asarray(values_half, dtype=float)
    N = vals.size
    thetas = np.linspace(0.0, np.pi, N, endpoint=False)
    thetas_pad = np.concatenate([thetas, [np.pi]])
    vals_pad   = np.concatenate([vals,   [vals[-1]]])

    f_half = interp1d(
        thetas_pad, vals_pad,
        kind=('linear' if kind == 'linear' else 'cubic'),
        assume_sorted=True, bounds_error=False, fill_value='extrapolate'
    )
    theta_full = np.linspace(0.0, 2.0*np.pi, n_full, endpoint=False)
    phi = theta_full % (2.0*np.pi)
    arg = np.where(phi < np.pi, phi, np.pi - (phi - np.pi)) 
    vals_full = f_half(arg)
    return theta_full, vals_full

def plot_engine_3d(
    designer,
    nx=500, ny=360,
    az_interp='cubic',
    axial_temp_interp='linear',
    log_temp=False,
    elev=25, azim=-60
):
    props = designer.engine.engineProps
    r_raw = np.asarray(props[:, 0], dtype=float)  
    x_raw = np.asarray(props[:, 1], dtype=float) 
    r_of_x = interp1d(x_raw, r_raw, kind='cubic', fill_value='extrapolate')

    x_min, x_max = float(x_raw.min()), float(x_raw.max())
    x_grid = np.linspace(x_min, x_max, nx)
    r_grid = r_of_x(x_grid)

    theta_full = np.linspace(0.0, 2.0*np.pi, ny, endpoint=False)
    X = np.repeat(x_grid[:, None], ny, axis=1)
    R = np.repeat(r_grid[:, None], ny, axis=1)
    Y = R * np.cos(theta_full)[None, :]
    Z = R * np.sin(theta_full)[None, :]

    nseg = int(designer.n_segments)

    full_rings = []
    for i in range(nseg):
        T_half = np.array(designer.T_wg_array[i, :], dtype=float)
        if np.all(~np.isfinite(T_half)):
            T_half = np.zeros_like(T_half)
        _, T_full = expand_half_to_full(T_half, n_full=ny, kind=az_interp)
        full_rings.append(T_full) 

    seg_x_bounds = np.array([get_segment_bounds(designer, i) for i in range(nseg)])
    seg_centers = seg_x_bounds.mean(axis=1)

    T_overlay = np.zeros((nx, ny), dtype=float)
    if axial_temp_interp == 'step':
        for i in range(nseg):
            x0, x1 = seg_x_bounds[i]
            mask = (x_grid >= x0) & (x_grid <= x1 if i == nseg - 1 else x_grid < x1)
            if not np.any(mask):
                continue
            Trow = full_rings[i]
            T_overlay[mask, :] = Trow[None, :]
    else:
        centers = seg_centers
        for ix, xval in enumerate(x_grid):
            idx = np.argsort(np.abs(centers - xval))
            i0, i1 = idx[0], (idx[1] if len(idx) > 1 else idx[0])
            c0, c1 = centers[i0], centers[i1]
            if i0 == i1 or c0 == c1:
                w0, w1 = 1.0, 0.0
            else:
                d0, d1 = abs(xval - c0), abs(xval - c1)
                w0 = d1 / (d0 + d1)
                w1 = d0 / (d0 + d1)
            T_overlay[ix, :] = w0 * full_rings[i0] + w1 * full_rings[i1]

    eps = 1e-9 
    temps = T_overlay.copy()
    if log_temp:
        temps = np.log10(np.maximum(temps, eps))

    tmin, tmax = np.nanmin(temps), np.nanmax(temps)
    if not np.isfinite(tmin) or not np.isfinite(tmax):
        raise RuntimeError("plot_engine_3d: temperature overlay invalid (all NaN?).")

    norm = Normalize(vmin=tmin, vmax=tmax)
    cmap = plt.get_cmap('inferno')
    facecolors = cmap(norm(temps))[:-1, :-1, :]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                    facecolors=facecolors, linewidth=0, antialiased=False, shade=False)

    ax.view_init(elev=elev, azim=azim)
    ax.set_xlabel('Axial x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z [m]')
    ax.set_title('Engine Hot-Wall Temperature ({} scale)'.format('log10' if log_temp else 'linear'))

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.12)
    cbar.set_label('log10(T_hot-wall [K])' if log_temp else 'T_hot-wall [K]')

    try:
        xr = X.max() - X.min()
        yr = Y.max() - Y.min()
        zr = Z.max() - Z.min()
        maxr = max(xr, yr, zr)
        ax.set_box_aspect((xr/maxr, yr/maxr, zr/maxr))
    except Exception:
        pass

    plt.tight_layout()
    plt.show()

def plot_segment_profiles(
    designer,
    seg_index,
    az_oversample=8,
    az_interp='cubic',
    pressure_units='Pa',
    log_temp=False
):

    i = int(seg_index)
    naz_half = int(designer.n_azimuthal)


    avg_props, _, _ = designer.get_segment_avg_props(i)
    r_seg = float(avg_props[0])
    try:
        chosen = designer.result['segments'][i]
    except Exception:
        raise RuntimeError("Run optimize() first to populate designer.result['segments'].")

    g = chosen['geometry']
    h = float(g.height)
    r_channel = r_seg + float(designer.inner_wall_thickness) + 0.5 * h

    n_fine = naz_half * az_oversample
    theta_half_fine = np.linspace(0.0, np.pi, n_fine, endpoint=False)

    def per_half(vals):
        return half_ring_interp(vals, kind=az_interp)(theta_half_fine)

    T_sat = per_half(designer.Tsat_array[i, :])
    T_c   = per_half(designer.T_c_array[i, :])
    T_wg  = per_half(designer.T_wg_array[i, :])
    T_fb  = per_half(designer.T_fb_array[i, :])
    P     = per_half(designer.P_array[i, :])


    L_total_half = np.pi * r_channel
    s = theta_half_fine / np.pi * L_total_half

    eps = 1e-9
    if log_temp:
        T_sat = np.log10(np.maximum(T_sat, eps))
        T_c   = np.log10(np.maximum(T_c,   eps))
        T_wg  = np.log10(np.maximum(T_wg,  eps))
        T_fb  = np.log10(np.maximum(T_fb,  eps))
        temp_ylabel = 'log10(Temperature [K])'
    else:
        temp_ylabel = 'Temperature [K]'

    if str(pressure_units).lower() == 'psi':
        P = P / 6894.757293 
        p_ylabel = 'Pressure [psi]'
    else:
        p_ylabel = 'Pressure [Pa]'


    fig, ax1 = plt.subplots()
    ax1.plot(s, T_sat, label='T_sat')
    ax1.plot(s, T_c,   label='T_coolant')
    ax1.plot(s, T_wg,  label='T_hot-wall')
    ax1.plot(s, T_fb,  label='T_fin-base')
    ax1.set_xlabel('Distance along channel s [m]')
    ax1.set_ylabel(temp_ylabel)

    ax2 = ax1.twinx()
    ax2.plot(s, P, linestyle='--', label='Pressure')
    ax2.set_ylabel(p_ylabel)


    l1, lab1 = ax1.get_legend_handles_labels()
    l2, lab2 = ax2.get_legend_handles_labels()
    ax1.legend(l1 + l2, lab1 + lab2, loc='best')

    ax1.grid(True, alpha=0.3)
    ax1.set_title(f'Segment {i}: Thermal & Pressure vs Channel Distance (Half Turn)')
    plt.tight_layout()
    plt.show()

