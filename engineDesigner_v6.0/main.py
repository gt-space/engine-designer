from EngineParameters import *
from Optimization.ContourCode import runContourCode
from Geometry import Contour, CoolingChannel
from Solvers import HotGasFlowSolver, ChannelSolver
from Outputs import Geometry, HotGasFlow, ChannelPerformance

import matplotlib.pyplot as plt


def main():

    # --------------------------------
    # Run contour optimization
    # --------------------------------
    runContourCode()

    # --------------------------------
    # Get engine contour
    # --------------------------------
    N = 120

    contour = Contour.getCurrentContour(N)

    # --------------------------------
    # Get hot gas flow
    # --------------------------------
    hot_gas = HotGasFlowSolver.getHotGasFlow(contour)

    # --------------------------------
    # Selected cooling channel geometry
    # --------------------------------
    N_ch = 90

    w_t = 1.05e-3       # throat channel width [m]
    AR_t = 2.6           # throat aspect ratio

    C_chamber = 2.0      # chamber width / throat width
    C_exit = 1.6         # exit width / throat width

    channel = CoolingChannel.getCoolingChannelGeometryFromPoints(
        contour,
        N_ch,
        w_t,
        AR_t,
        C_chamber,
        C_exit
    )

    # --------------------------------
    # Solve coolant flow
    # --------------------------------
    coolant_flow = ChannelSolver.CoolantChannelSolver(
        T_cool_in,
        contour,
        channel,
        hot_gas
    )

    # --------------------------------
    # Print key performance values
    # --------------------------------
    dP = coolant_flow["P_cool"][-1] - coolant_flow["P_cool"][0]
    dP_psi = dP * 0.000145038

    max_Twg = max(coolant_flow["Twg"])

    throat_index = contour["throat_index"]

    print("\nSelected Cooling Channel Geometry")
    print("--------------------------------")
    print(f"Number of channels:       {N_ch}")
    print(f"Throat channel width:     {w_t*1000:.3f} mm")
    print(f"Throat aspect ratio:      {AR_t:.2f}")
    print(f"Throat channel height:    {channel['h_ch'][throat_index]*1000:.3f} mm")
    print(f"Throat fin width:         {channel['w_fin'][throat_index]*1000:.3f} mm")
    print(f"C_chamber:                {C_chamber:.2f}")
    print(f"C_exit:                   {C_exit:.2f}")

    print("\nPerformance")
    print("--------------------------------")
    print(f"Coolant pressure drop:    {dP_psi:.2f} psi")
    print(f"Maximum Twg:              {max_Twg:.1f} K")
    print(f"Minimum fin width:        {min(channel['w_fin'])*1000:.3f} mm")

    # --------------------------------
    # Outputs
    # --------------------------------
    Geometry.plotGeometry(
        contour,
        channel
    )

    HotGasFlow.plotHotGasFlow(
        contour,
        hot_gas
    )

    ChannelPerformance.plotChannelPerformance(
        contour,
        hot_gas,
        coolant_flow
    )

    # Show all figures together
    plt.show()


if __name__ == "__main__":
    main()