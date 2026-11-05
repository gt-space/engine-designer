from EngineParameters import *
from Optimization.ContourOptimization import setContourWeights, setContourBounds, optimizeContour
from Optimization import ManifoldOptimization
from Geometry import Contour, CoolingChannel, Manifold
from Solvers import HotGasFlowSolver, ChannelSolver
from Outputs import Geometry, HotGasFlow, ChannelPerformance, ManifoldPerformance

import matplotlib.pyplot as plt


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

    # 5. Solve coolant flow along channels
    coolant_flow = ChannelSolver.CoolantChannelSolver(
        T_cool_in, contour, channel, hot_gas
    )

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

    print("\nSelected Cooling Channel & Manifold Summary")
    print("-------------------------------------------")
    print(f"Number of channels:          {N_ch}")
    print(f"Throat channel width:        {w_t*1000:.3f} mm")
    print(f"Throat channel height:       {channel['h_ch'][throat_index]*1000:.3f} mm")
    print(f"Channel pressure drop:       {dP_psi:.2f} psi")
    print(f"Maximum Twg:                 {max_Twg:.1f} K")
    print(f"Manifold Inlet Split Angle:  {manifold_opt['inlet_angle_deg']:.2f} deg")
    print(f"Manifold Inlet Inner Radius: {manifold_opt['r_manifold_profile'][0]*1e3:.2f} mm")
    print(f"Manifold 180° End Radius:    {manifold_opt['r_manifold_profile'][-1]*1e3:.2f} mm")
    print(f"Channel Flow Maldistribution:{manifold_opt['maldistribution_pct']:.3f} %")
    print(f"Total Feed In Pressure Req:  {manifold_opt['P_inlet_required_Pa']*0.000145038:.2f} psi")

    # 8. Plot all outputs
    Geometry.plotGeometry(contour, channel)
    HotGasFlow.plotHotGasFlow(contour, hot_gas)
    ChannelPerformance.plotChannelPerformance(contour, hot_gas, coolant_flow)
    ManifoldPerformance.plotManifoldResults(contour, manifold_opt, cad_loft)

    plt.show()


if __name__ == "__main__":
    main()