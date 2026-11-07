import numpy as np
from pathlib import Path
from EngineParameters import *
from Optimization.ContourOptimization import setContourWeights, setContourBounds, optimizeContour
from Optimization import ManifoldOptimization, LinerShellOptimization
from Geometry import Contour, CoolingChannel, Manifold, LinerShell
from Solvers import HotGasFlowSolver, ChannelSolver, StructuralSolver
from Outputs import Geometry as PlotGeometry, HotGasFlow, ChannelPerformance, ManifoldPerformance, Structural

import matplotlib.pyplot as plt

#Exports 3-column (X, Y, Z in mm) headerless curve files directly compatible with SolidWorks
def exportSolidWorksCurveFiles(contour, liner_geom, shell_geom):
    geometry_dir = Path(__file__).resolve().parent / "Geometry"
    geometry_dir.mkdir(parents=True, exist_ok=True)

    z_mm = contour["z"] * 1e3
    zeros = np.zeros_like(z_mm)

    curves = {
        geometry_dir / "SW_Curve_1_HotGas_Inner.txt": contour["r"] * 1e3,
        geometry_dir / "SW_Curve_2_Channel_Floor.txt": liner_geom["r_liner_outer"] * 1e3,
        geometry_dir / "SW_Curve_3_Channel_Roof_Shell_Inner.txt": shell_geom["r_shell_inner"] * 1e3,
        geometry_dir / "SW_Curve_4_Outer_Shell_Wall.txt": shell_geom["r_shell_outer"] * 1e3,
    }

    for filepath, r_mm in curves.items():
        xyz_pts = np.column_stack((r_mm, zeros, z_mm))
        np.savetxt(str(filepath), xyz_pts, fmt="%.5f", delimiter="\t")
    print(f"[Geometry] Saved 4 SolidWorks XYZ Curve files to {geometry_dir}")


def main():
    # 1. Run contour optimization
    weights = setContourWeights()
    bounds = setContourBounds()
    optimizeContour(weights, bounds)

    # 2. Get engine contour
    N = 120
    contour = Contour.getCurrentContour(N)

    # 3. Get hot gas flow
    hot_gas = HotGasFlowSolver.getHotGasFlow(contour)

    # 4. Selected cooling channel geometry
    N_ch = 90
    w_t = 1.05e-3       # throat channel width [m]
    AR_t = 2.6          # throat aspect ratio
    C_chamber = 2.0     # chamber width / throat width
    C_exit = 1.6        # exit width / throat width

    channel = CoolingChannel.getCoolingChannelGeometryFromPoints(
        contour, N_ch, w_t, AR_t, C_chamber, C_exit
    )

    # 5. Solve initial coolant flow + coupled Smooth Liner & Outer Shell Thicknesses
    coolant_flow_init = ChannelSolver.CoolantChannelSolver(
        T_cool_in, contour, channel, hot_gas
    )

    struct_results = LinerShellOptimization.optimizeLinerAndShellProfiles(
        contour, channel, hot_gas, coolant_flow_init,
        liner_material=liner_material_name,
        shell_material=shell_material_name
    )

    # Use updated channel (with variable priority liner thickness tw) and updated coolant flow
    channel = struct_results["channel_updated"]
    coolant_flow = struct_results["coolant_flow_updated"]

    # Export Full Liner + Channel + Outer Shell CAD Table & SolidWorks XYZ Curves
    LinerShell.exportLinerShellCADToCSV(
        contour, channel, struct_results["liner"], struct_results["shell"]
    )
    exportSolidWorksCurveFiles(contour, struct_results["liner"], struct_results["shell"])

    # 6. Optimize Single-Inlet Tapered Manifold coupled to ChannelSolver results
    manifold_channel_geom = Manifold.extractManifoldChannelInterface(contour, channel, coolant_flow)
    manifold_struct_params = Manifold.setManifoldSizeLimitsFromContour(contour, channel)
    manifold_engine_params = Manifold.loadEngineCoolantParameters(coolant_flow)

    manifold_opt = ManifoldOptimization.optimizeManifold(
        engine_params=manifold_engine_params,
        channel_geom=manifold_channel_geom,
        struct_params=manifold_struct_params,
        n_ctrl=5,
        shape="circular",
        max_iter=50
    )

    cad_loft = Manifold.build3DCADManifoldLoft(manifold_opt, manifold_channel_geom)
    Manifold.exportManifoldGeometryToCSV(cad_loft, manifold_opt)

    # 7. Print key performance values
    dP = coolant_flow["P_cool"][-1] - coolant_flow["P_cool"][0]
    dP_psi = dP * 0.000145038
    max_Twg = max(coolant_flow["Twg"])
    throat_index = contour["throat_index"]

    print("\nSelected Cooling Channel, Liner, Shell & Manifold Summary")
    print("---------------------------------------------------------")
    print(f"Number of channels:          {N_ch}")
    print(f"Throat channel width:        {w_t*1000:.3f} mm")
    print(f"Throat channel height:       {channel['h_ch'][throat_index]*1000:.3f} mm")
    print(f"Liner thickness (throat):    {channel['tw'][throat_index]*1000:.3f} mm ({struct_results['liner']['governing_mode'][throat_index]})")
    print(f"Liner thickness (max):       {max(channel['tw'])*1000:.3f} mm")
    print(f"Outer Shell (throat):        {struct_results['shell']['t_shell_selected'][throat_index]*1000:.3f} mm ({struct_results['shell']['governing_mode'][throat_index]})")
    print(f"Outer Shell (max):           {max(struct_results['shell']['t_shell_selected'])*1000:.3f} mm")
    print(f"Total Chamber Metal Mass:    {struct_results['mass_total_chamber_kg']:.2f} kg")
    print(f"Channel pressure drop:       {dP_psi:.2f} psi")
    print(f"Maximum Twg:                 {max_Twg:.1f} K")
    print(f"Manifold Inlet Split Angle:  {manifold_opt['inlet_angle_deg']:.2f} deg")
    print(f"Manifold Inlet Inner Radius: {manifold_opt['r_manifold_profile'][0]*1e3:.2f} mm")
    print(f"Manifold 180° End Radius:    {manifold_opt['r_manifold_profile'][-1]*1e3:.2f} mm")
    print(f"Channel Flow Maldistribution:{manifold_opt['maldistribution_pct']:.3f} %")
    print(f"Total Feed In Pressure Req:  {manifold_opt['P_inlet_required_Pa']*0.000145038:.2f} psi")

    # 8. Plot all outputs
    PlotGeometry.plotGeometry(contour, channel)
    HotGasFlow.plotHotGasFlow(contour, hot_gas)
    ChannelPerformance.plotChannelPerformance(contour, hot_gas, coolant_flow)
    Structural.plotStructuralSizing(contour, struct_results)
    ManifoldPerformance.plotManifoldResults(contour, manifold_opt, cad_loft)

    plt.show()


if __name__ == "__main__":
    main()