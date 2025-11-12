import CoolProp.CoolProp as CP
import numpy as np
import pandas as pd
from bartz import bartz
from typing import Tuple, Dict, List
from scipy.optimize import root
from scipy.interpolate import interp1d
from dataclasses import dataclass

# Data classes and structures
print("LOADED calChamber2.py")

@dataclass
class Geometry:
    height: float
    width: float
    num_channels: int

class MaterialProperties:
    def __init__(self, conductivity: float, yield_strength: float, ultimate_strength: float):
        self.k = conductivity
        self.yield_strength = yield_strength
        self.ultimate_strength = ultimate_strength

class CoolingChannelDesigner:
    def __init__(self, engine, material: 'MaterialProperties', settings: Dict):
        self.engine = engine
        self.material = material
        self.settings = settings
        self.n_segments = settings['n_segments']
        self.n_azimuthal = settings['n_azimuthal_divisions']
        self.P_inlet = settings['P_inlet']  # in Pascals
        self.T_inlet = settings['T_in']
        self.wall_temp_profile = np.full(self.n_segments, settings['Target_Wall_Temp'])
        self.inner_wall_thickness = settings['inner_wall_thickness']
        self.result = None
        self.T_wg_array = np.full((self.n_segments, self.n_azimuthal), np.nan)
        self.T_fb_array = np.full((self.n_segments, self.n_azimuthal), np.nan)
        self.T_c_array = np.full((self.n_segments, self.n_azimuthal), np.nan)
        self.q_array = np.full((self.n_segments, self.n_azimuthal), np.nan)
        self.P_array = np.full((self.n_segments, self.n_azimuthal), np.nan)
        self.mdot_array = np.full((self.n_segments, self.n_azimuthal), np.nan)
        self.superheat = settings['Superheat']
        self.Tsat_array = np.full((self.n_segments, self.n_azimuthal), np.nan)

    def optimize(self):
        converged = False
        iteration = 0

        segment_candidates = self.axial_solver(
            T_wall_target_profile=self.wall_temp_profile,
            geometry_controls=self.settings['geometry_controls'],
            n_segments=self.n_segments,
            n_azimuthal=self.settings['n_azimuthal_divisions'],
            sweep_ranges=self.settings['sweep_ranges'],
            weights=self.settings['weights'],
        )

        self.segment_candidates = segment_candidates

        while not converged:
            segment_results = []
            for i in range(self.n_segments):
                candidates = segment_candidates[i]
                for candidate in candidates:
                    result = {'segments': segment_results + [candidate] + [segment_candidates[j][0] for j in range(i+1, self.n_segments)]}
                    updated_profile, FOS_ok = self.structural_solver(result, self.material, self.settings)
                    if FOS_ok:
                        segment_results.append(candidate)
                        break
            self.result = {'segments': segment_results}
            updated_profile, FOS_ok = self.structural_solver(self.result, self.material, self.settings)
            converged = FOS_ok and np.allclose(self.wall_temp_profile, updated_profile, rtol=0.01)
            self.wall_temp_profile = updated_profile
            iteration += 1

        return {
            'segments': segment_results,
            'candidates': segment_candidates
        }

    def summarize(self):
        if not self.result or not self.result.get('segments'):
            print("No results available.")
            return

        x = np.asarray(self.engine.engineProps[:, 1], dtype=float)
        total_engine_length = float(np.max(x) - np.min(x))
        L_ax = total_engine_length / float(self.n_segments)

        total_length = 0.0
        total_channels = 0
        total_mdot = 0.0

        print("Calorimetry Summary")
        for seg in self.result['segments']:
            i = int(seg['segment'])
            g = seg['geometry']

            fin_width = (L_ax - (g.width * g.num_channels)) / max(1, g.num_channels)

            seg_total_mdot = seg['mdot'] * 2.0

            total_length += L_ax
            total_channels += g.num_channels
            total_mdot += seg_total_mdot

            print(
                f"Segment {i}: "
                f"Fin width={fin_width*39.37:.6f} in, "
                f"Fin height={g.height*39.37:.6f} in, "
                f"Channel width={g.width*39.37:.6f} in, "
                f"Number of Channels={g.num_channels:d}, "
            )

        print("\nTotals")
        print(f"Total engine length = {total_length*39.37:.6f} in")
        print(f"Total number of channels = {total_channels:d}")
        print(f"Total mdot = {total_mdot:.6f} kg/s")

    def export_results(self, filename="engine_sketch_points.csv", nx=1500,
                            min_points_per_fin_top=8, min_points_per_side=2):
        import csv

        props = self.engine.engineProps
        x_raw = np.asarray(props[:, 1], dtype=float)
        r_raw = np.asarray(props[:, 0], dtype=float)


        r_of_x = interp1d(x_raw, r_raw, kind='cubic', fill_value='extrapolate')
        x_grid = np.linspace(x_raw.min(), x_raw.max(), int(nx))
        r_grid = r_of_x(x_grid)
        drdx = np.gradient(r_grid, x_grid)


        ds = np.sqrt(1.0 + drdx**2) * np.gradient(x_grid)
        s_grid = np.concatenate([[0.0], np.cumsum(ds[1:])])

        inv_norm = np.sqrt(1.0 + drdx**2)
        nx_hat = -drdx / inv_norm
        nr_hat =  1.0  / inv_norm

        r_on = interp1d(x_grid, r_grid,  kind='linear', assume_sorted=True, fill_value='extrapolate')
        nx_on = interp1d(x_grid, nx_hat, kind='linear', assume_sorted=True, fill_value='extrapolate')
        nr_on = interp1d(x_grid, nr_hat, kind='linear', assume_sorted=True, fill_value='extrapolate')

        def offset_curve(x_pts, offset):
            nxp = nx_on(x_pts); nrp = nr_on(x_pts)
            return x_pts + offset * nxp, r_on(x_pts) + offset * nrp

        rows = []

        for idx, (xv, rv) in enumerate(zip(x_grid, r_grid)):
            rows.append(["contour", -1, -1, idx, float(xv), float(rv)])

        x_in, r_in = offset_curve(x_grid, float(self.inner_wall_thickness))
        for idx, (xv, rv) in enumerate(zip(x_in, r_in)):
            rows.append(["inner_wall", -1, -1, idx, float(xv), float(rv)])

    # Objective function
    def objective_function(self, dP: float, delta_T: float, mdot: float,
                           T_wall: float, Tw_target: float, weights: Dict, sweep_ranges) -> float:
        return (
            weights['dp'] * dP / self.P_inlet +
            weights['dT'] * delta_T +
            weights['mdot'] * mdot / np.max(sweep_ranges['mdot']) +
            weights['Tw'] * (T_wall / Tw_target)
        )

    # Geometric constraint checker
    def violates_constraints(self, h: float, w: float, N: int, mdot: float, fin_width: float, aspect_ratio: float, geom_ctrl: Dict) -> bool:
        fix = geom_ctrl.get('fix', [])
        constraints = geom_ctrl.get('constraints', {})

        if 'h' in fix and not np.isclose(h, geom_ctrl['h']):
            return True
        if 'w' in fix and not np.isclose(w, geom_ctrl['w']):
            return True
        if 'N' in fix and not N == geom_ctrl['N']:
            return True
        if 'mdot' in fix and not np.isclose(mdot, geom_ctrl['mdot']):
            return True

        if 'min_fin_thickness' in constraints and fin_width < constraints['min_fin_thickness']:
            return True
        if 'max_aspect_ratio' in constraints and aspect_ratio > constraints['max_aspect_ratio']:
            return True
        if 'min_h' in constraints and h < constraints['min_h']:
            return True
        if 'max_h' in constraints and h > constraints['max_h']:
            return True
        if 'min_w' in constraints and w < constraints['min_w']:
            return True
        if 'max_w' in constraints and w > constraints['max_w']:
            return True

        return False

    # Axial solver
    def axial_solver(self, T_wall_target_profile, geometry_controls, n_segments, sweep_ranges, weights, n_azimuthal):
        segment_candidates = {}
        segment_results = []

        for i in range(n_segments):
            candidates = []
            Tw_target = T_wall_target_profile[i]
            geom_ctrl = geometry_controls[i]
            avg_props, start, end = self.get_segment_avg_props(i)
            r = avg_props[0]
            L_ax = self.engine.engineProps[199, 1] / n_segments
            A_wg = (np.pi * r / n_azimuthal) * L_ax

            n_combinations_tried = 0  # DEBUG
            n_combinations_valid = 0  # DEBUG

            for h in sweep_ranges['h']:
                for w in sweep_ranges['w']:
                    for N in sweep_ranges['N']:
                        for mdot in sweep_ranges['mdot']:

                            D_h = 2 * h * w / (h + w)
                            fin_width = (L_ax - (w * N)) / N
                            aspect_ratio = h / fin_width if fin_width > 0 else np.inf
                            r_channel = r + self.inner_wall_thickness + h / 2
                            L_az_channel = (np.pi * r_channel) / n_azimuthal
                            A_flow = N * w * h
                            dV = A_flow * L_az_channel

                            n_combinations_tried += 1  # DEBUG

                            if self.violates_constraints(h, w, N, mdot, fin_width, aspect_ratio, geom_ctrl):
                                continue

                            Tc_local = self.T_inlet
                            Pc_local = self.P_inlet
                            dT_total = 0
                            dP_total = 0

                            terminate_geometry = False

                            T_wg_list = []
                            T_fb_list = []
                            Tc_list = []
                            q_list = []
                            Pc_list = []
                            mdot_list = []
                            Tsat_list = []

                            for j in range(n_azimuthal):
                                try:
                                    rho = CP.PropsSI('D', 'T', Tc_local, 'P', Pc_local, "Water")
                                except ValueError:
                                    terminate_geometry = True
                                    break

                                mu = CP.PropsSI('V', 'T', Tc_local, 'P', Pc_local, "Water") / rho
                                Cp = CP.PropsSI('C', 'T', Tc_local, 'P', Pc_local, "Water")
                                cond_c = CP.PropsSI('CONDUCTIVITY', 'T', Tc_local, 'P', Pc_local, "Water")
                                v = mdot / (rho * A_flow)
                                dM = rho * dV
                                dt = L_az_channel / v

                                Re = (v * D_h) / mu
                                Pr = mu * rho * Cp / cond_c
                                f = self.get_friction(Re, D_h)

                                result = self.solve_wall_temperature(i, Tc_local, h, w, N, A_wg, L_az_channel, fin_width,
                                                                     f, Re, Pr, cond_c, D_h, mu, Tw_target)

                                if result is None:
                                    terminate_geometry = True
                                    break
                                T_wall, T_fin_base, Q_in = result

                                dP = self.pressureDrop(rho, L_az_channel, D_h, v, f)
                                Pc_local -= dP
                                dP_total += dP

                                dT = (Q_in * dt) / (dM * Cp)
                                Tc_local += dT
                                dT_total += dT

                                try:
                                    Tsat_local = self.get_Tsat_from_pressure(Pc_local)
                                except ValueError:
                                    terminate_geometry = True
                                    break
                                
                                superheat = T_fin_base - Tsat_local
                                if superheat > self.superheat:
                                    terminate_geometry = True
                                    break

                                T_wg_list.append(T_wall)
                                T_fb_list.append(T_fin_base)
                                Tc_list.append(Tc_local)
                                q_list.append(Q_in)
                                Pc_list.append(Pc_local)
                                mdot_list.append(mdot)
                                Tsat_list.append(Tsat_local)

                            if terminate_geometry:
                                continue

                            max_dT = self.settings.get('max_dT')
                            max_dP = self.settings.get('max_dP')
                            if (max_dT is not None and dT_total > max_dT) or \
                               (max_dP is not None and dP_total > max_dP):
                                continue

                            n_combinations_valid += 1  # DEBUG

                            score = self.objective_function(dP_total, dT_total, mdot, T_wall, Tw_target, weights, sweep_ranges)

                            candidates.append({
                                'segment': i,
                                'Tw_target': Tw_target,
                                'Tc': Tc_local,
                                'dP': dP_total,
                                'geometry': Geometry(height=h, width=w, num_channels=N),
                                'mdot': mdot,
                                'score': score,
                                'dT': dT_total,
                                'T_wall_exit': T_wall,
                                'T_wg_list': T_wg_list,
                                'T_fb_list': T_fb_list,
                                'T_c_list':  Tc_list,
                                'q_list':    q_list,
                                'P_list':    Pc_list,
                                'mdot_list': mdot_list,
                                'Tsat_list': Tsat_list,
                            })

            if candidates:
                best = min(candidates, key=lambda c: (c['score']))
                segment_candidates[i] = [best]
            else:
                print(f"No valid configurations for segment {i}, skipping this segment.")
                segment_candidates[i] = []
                continue

            self.T_wg_array[i, :] = np.asarray(best['T_wg_list'], dtype=float)
            self.T_fb_array[i, :] = np.asarray(best['T_fb_list'], dtype=float)
            self.T_c_array[i, :]  = np.asarray(best['T_c_list'],  dtype=float)
            self.q_array[i, :]    = np.asarray(best['q_list'],     dtype=float)
            self.P_array[i, :]    = np.asarray(best['P_list'],     dtype=float)
            self.mdot_array[i, :] = np.asarray(best['mdot_list'],  dtype=float)
            self.Tsat_array[i, :] = np.asarray(best['Tsat_list'],  dtype=float)

            g = best['geometry']
            print(f"Solved segment {i + 1} — mdot: {best['mdot']:.4f} kg/s, T_wall_exit: {best['T_wall_exit']:.2f} K, "
                  f"dP: {best['dP']:.2f} Pa, geometry: h={g.height:.4f}, w={g.width:.4f}, N={g.num_channels}, "
                  f"Tc_exit={best['Tc']:.2f}, dT={best['dT']:.2f}")
            print(f"Segment {i} — {n_combinations_tried} combinations tried, {n_combinations_valid} valid")

        unsolved_segments = [i for i, cands in segment_candidates.items() if not cands]
        if unsolved_segments:
            raise RuntimeError(f"Axial solver failed for segment(s): {unsolved_segments}. "
                               f"No valid configurations found, check constraints.")
        return segment_candidates

    # Structural Solver
    def structural_solver(self, result: Dict, material: 'MaterialProperties', settings: Dict) -> Tuple[np.ndarray, bool]:
        T_wall_target_profile = np.array([seg['Tw_target'] for seg in result['segments']])
        coolant_temp_profile = np.array([seg['Tc'] for seg in result['segments']])

        updated_T_wall_target = []
        FOS_check = True

        for Tw_target, Tc in zip(T_wall_target_profile, coolant_temp_profile):
            stress = np.nan
            # FOS = material.yield_strength / stress
            FOS = settings['min_FOS']

            if FOS < settings['min_FOS']:
                FOS_check = False
                Tw_target_new = Tw_target * 0.95
            else:
                Tw_target_new = Tw_target

            updated_T_wall_target.append(Tw_target_new)

        return np.array(updated_T_wall_target), FOS_check

    # Heat transfer stuff
    def solve_wall_temperature(self, i, Tc_local, h, w, N, A_wg, L_az_channel, fin_width,
                               f, Re, Pr, cond_c, D_h, mu, T_max_allowed) -> Tuple[float, float, float] | None:

        avg_props, start, end = self.get_segment_avg_props(i)
        R1 = avg_props[0]

        k = self.material.k
        T_cb = Tc_local  # Local coolant bulk temp
        h_c = self.Gnielinski(f, Re, Pr, cond_c, D_h, mu)

        m = np.sqrt((2 * h_c * fin_width) / k)
        n_f = (np.tanh(m * (h / fin_width))) / (m * (h / fin_width))
        h_cf = h_c * ((w + 2 * n_f * h) / (w + fin_width))

        def residual(T_wg):
            h_g, q_in, T_aw = self.average_bartz(self.engine, T_wg, start, end)
            Q_conv = h_g * A_wg * (T_aw - T_wg)
            T_fb = T_wg - (Q_conv * (self.inner_wall_thickness / (A_wg * k)))
            Q_cool = h_cf * (L_az_channel * N * (2 * h + w)) * (T_fb - T_cb)
            return Q_conv - Q_cool

        sol = root(residual, T_max_allowed, method='hybr')
        if not sol.success:
            return None

        T_wg = float(sol.x[0])

        if T_wg > T_max_allowed:
            return None

        h_g, q_in, T_aw = self.average_bartz(self.engine, T_wg, start, end)
        Q_conv = h_g * A_wg * (T_aw - T_wg)
        T_fb = T_wg - (Q_conv * (self.inner_wall_thickness / (A_wg * k)))
        return T_wg, T_fb, Q_conv

    # Other functions
    def pressureDrop(self, rho, L, D_h, v, f) -> float:
        return f * (L / D_h) * 0.5 * rho * v**2

    def get_segment_avg_props(self, i: int) -> np.ndarray:
        n_props = self.engine.engineProps.shape[0]
        start = int(i * n_props / self.n_segments)
        end = int((i + 1) * n_props / self.n_segments)
        return np.mean(self.engine.engineProps[start:end, :], axis=0), start, end

    def average_bartz(self, engine, T_wg: float, start: int, end: int) -> Tuple[float, float, float]:
        # Averages Bartz parameters over a segment of engine stations.
        h_g_list = []
        q_in_list = []
        T_aw_list = []

        for station in range(start, end):
            h_g, q_in, T_aw = bartz(engine, T_wg, station)
            h_g_list.append(h_g)
            q_in_list.append(q_in)
            T_aw_list.append(T_aw)

        return (
            np.mean(h_g_list),
            np.mean(q_in_list),
            np.mean(T_aw_list)
        )

    def get_friction(self, Re, D_hyd):
        def colebrook(Re, D_hyd):
            k = 0.005e-3  # Surface roughness in meters (5 microns)
            fric = 1 / ((-1.8 * np.log10((k / D_hyd / 3.7) ** 1.11 + 6.9 / Re)) ** 2)
            cole_diff = -2 * np.log10(k / D_hyd / 3.7 + 2.51 / Re / np.sqrt(fric)) - 1 / np.sqrt(fric)
            while np.abs(cole_diff) > 1e-5:
                deriv = -2 / (np.log(10) * (k / D_hyd / 3.7 + 2.51 / Re / np.sqrt(fric))) * \
                        (2.51 / Re * (-0.5 * fric ** (-1.5))) + 0.5 * fric ** (-1.5)
                fric = fric - cole_diff / deriv
                cole_diff = -2 * np.log10(k / D_hyd / 3.7 + 2.51 / Re / np.sqrt(fric)) - 1 / np.sqrt(fric)
            return fric

        return 64 / Re if Re <= 2320 else colebrook(Re, D_hyd)

    def Gnielinski(self, f, Re, Pr, cond_c, Dh, mu):
        Nu = (f / 8) * (Re - 1000) * Pr / (1 + 12.7 * ((f / 8) ** 0.5) * (Pr ** (2 / 3) - 1))
        hc = Nu * cond_c / Dh * (1000000 * mu / 0.2) ** 0.11
        return hc

    def get_Tsat_from_pressure(self, P: float) -> float:
        return float(CP.PropsSI('T', 'P', P, 'Q', 0.0, 'Water'))
