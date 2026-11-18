import numpy as np
from fullflow import Network, State, Component, Balance, SteadyState, DarcyWeisbach
from Functions import PressureDrop as pd

#Single feed pipe splitting into two symmetric 180-deg manifold streams at split angle alpha
class ManifoldInletSplit(Component):
    def __init__(self, name, network, mdot_pipe, P_pipe, P_man_inlet,
                 rho, A_pipe, A_man_inlet, split_angle_deg=90.0):
        self.name = name
        self.network = network
        self.mdot_pipe = mdot_pipe
        self.P_pipe = P_pipe
        self.P_man_inlet = P_man_inlet
        self.rho = rho
        self.A_pipe = A_pipe
        self.A_man_inlet = A_man_inlet
        self.split_angle_deg = split_angle_deg
        self.residual = 0.0
        self.setup()

    def evaluate_states(self):
        _, dP_static = pd.dPinletSplit(
            self.mdot_pipe.value,
            self.rho.value,
            self.A_pipe.value,
            self.A_man_inlet.value,
            self.split_angle_deg.value
        )
        self.residual = ((self.P_pipe.value - dP_static) - self.P_man_inlet.value) / 1e5

    @property
    def balances(self):
        return [(self.P_man_inlet, self.residual)]

#Circumferential manifold segment between channel i and channel i+1.
# Balances wall friction against deceleration momentum static pressure regain.
class ManifoldRingSegment(Component):
    def __init__(self, name, network, P_up, P_down, mdot_up, mdot_down,
                 A_up, A_down, Dh_up, Dh_down, dL, rho, mu, roughness):
        self.name = name
        self.network = network
        self.P_up = P_up
        self.P_down = P_down
        self.mdot_up = mdot_up
        self.mdot_down = mdot_down
        self.A_up = A_up
        self.A_down = A_down
        self.Dh_up = Dh_up
        self.Dh_down = Dh_down
        self.dL = dL
        self.rho = rho
        self.mu = mu
        self.roughness = roughness
        self.p_seg_residual = 0.0
        self.setup()

    def evaluate_states(self):
        dP_static_change, _ = pd.dPmanifoldSegment(
            self.mdot_up.value,
            self.mdot_down.value,
            self.A_up.value,
            self.A_down.value,
            self.Dh_up.value,
            self.Dh_down.value,
            self.dL.value,
            self.rho.value,
            self.mu.value,
            self.roughness.value
        )
        self.p_seg_residual = ((self.P_up.value + dP_static_change) - self.P_down.value) / 1e5

    @property
    def balances(self):
        return [(self.P_down, self.p_seg_residual)]

#90-degree branch takeoff + area contraction from manifold station i into cooling channel i.
class ManifoldBranchNode(Component):
    def __init__(self, name, network, P_man_local, P_ch_entry, mdot_man_in,
                 mdot_ch, A_man, A_ch, rho):
        self.name = name
        self.network = network
        self.P_man_local = P_man_local
        self.P_ch_entry = P_ch_entry
        self.mdot_man_in = mdot_man_in
        self.mdot_ch = mdot_ch
        self.A_man = A_man
        self.A_ch = A_ch
        self.rho = rho
        self.p_branch_residual = 0.0
        self.setup()

    def evaluate_states(self):
        dP_branch = pd.dPchannelBranch(
            self.mdot_man_in.value,
            self.mdot_ch.value,
            self.A_man.value,
            self.A_ch.value,
            self.rho.value
        )
        self.p_branch_residual = ((self.P_man_local.value - dP_branch) - self.P_ch_entry.value) / 1e5

    @property
    def balances(self):
        return [(self.P_ch_entry, self.p_branch_residual)]

# Builds and solves the 1D FullFlow network for a given manifold area schedule.
def solveManifoldFlowNetwork(A_man_profile, Dh_man_profile, inlet_angle_deg, engine_params, channel_geom):
    net = Network("Manifold_Cooling_Network")

    N_half = channel_geom["N_channels"] // 2
    mdot_total = engine_params["mdot_total"]
    mdot_half_target = 0.5 * mdot_total
    rho = engine_params["rho"]
    mu = engine_params["mu"]
    roughness = engine_params["roughness"]
    P_exit = engine_params["P_collector_exit"]

    # Dynamically compute channel friction factor from current (nominal or uncertain) parameters
    Dh_ch = channel_geom["Dh_channel"]
    scale_factor = channel_geom.get("scale_factor", 1.0)
    A_ch_single = channel_geom["A_channel_avg"]
    mdot_ch_single_nom = (mdot_half_target / N_half) / scale_factor

    Dh_nom = channel_geom.get("Dh_channel_nominal", Dh_ch)
    area_dh_scale = (Dh_ch / Dh_nom) ** 2
    A_ch_eff = A_ch_single * area_dh_scale
    A_ch_inlet_eff = channel_geom["A_ch_inlet"] * area_dh_scale

    V_ch_nom = mdot_ch_single_nom / (rho * (A_ch_single * area_dh_scale))
    Re_ch_nom = (rho * V_ch_nom * Dh_ch) / mu
    f_ch_base = float(pd.darcyWeishbachFriction(Re_ch_nom, roughness, Dh_ch))
    # Multiply by thermal/curvature multiplier from nominal coupled ChannelSolver if available
    f_mult = channel_geom.get("f_Save_multiplier", 1.0)
    f_channel_active = np.clip(f_ch_base * f_mult, 0.005, 0.35)

    d_theta = np.pi / N_half
    dL_segment = channel_geom["R_manifold"] * d_theta

    P_inlet = State(engine_params["P_inlet_guess"], bounds=(P_exit + 1e4, 5e8), keep_feasible=True)
    P_man_0 = State(engine_params["P_inlet_guess"] * 0.98, bounds=(P_exit + 1e3, 5e8), keep_feasible=True)

    ManifoldInletSplit(
        "Inlet_Split", net,
        mdot_pipe=mdot_total,
        P_pipe=P_inlet,
        P_man_inlet=P_man_0,
        rho=rho,
        A_pipe=engine_params["A_feed_pipe"],
        A_man_inlet=A_man_profile[0],
        split_angle_deg=inlet_angle_deg
    )

    P_nodes = [P_man_0]
    mdot_channels = []
    mdot_man_rem = State(mdot_half_target)

    for i in range(N_half):
        P_ch_in = State(P_man_0.value * 0.95, bounds=(P_exit + 100.0, 5e8), keep_feasible=True)
        # mdot_ch is now a direct assignable State representing ONE single channel's flow
        mdot_ch = State(mdot_ch_single_nom, bounds=(1e-6, mdot_half_target / scale_factor), keep_feasible=True)
        mdot_channels.append(mdot_ch)

        ManifoldBranchNode(
            f"Branch_{i}", net,
            P_man_local=P_nodes[i],
            P_ch_entry=P_ch_in,
            mdot_man_in=mdot_man_rem,
            mdot_ch=mdot_ch,
            A_man=A_man_profile[i],
            A_ch=A_ch_inlet_eff,
            rho=rho
        )

        DarcyWeisbach(
            f"CoolingChannel_{i}", net,
            mass_flow=mdot_ch,
            upstream_pressure=P_ch_in,
            downstream_pressure=P_exit,
            length=channel_geom["L_channel"],
            hydraulic_diameter=Dh_ch,
            density=rho,
            cross_sectional_area=A_ch_eff,
            friction_factor=f_channel_active
        )

        if i < N_half - 1:
            # Multiply single-channel flow by scale_factor for the lumped manifold bleed-off
            mdot_man_next = State.maximum(mdot_man_rem - mdot_ch * scale_factor, 1e-6)
            P_next = State(P_nodes[i].value, bounds=(P_exit + 1e3, 5e8), keep_feasible=True)
            P_nodes.append(P_next)

            ManifoldRingSegment(
                f"RingSeg_{i}", net,
                P_up=P_nodes[i],
                P_down=P_next,
                mdot_up=mdot_man_rem,
                mdot_down=mdot_man_next,
                A_up=A_man_profile[i],
                A_down=A_man_profile[i + 1],
                Dh_up=Dh_man_profile[i],
                Dh_down=Dh_man_profile[i + 1],
                dL=dL_segment,
                rho=rho,
                mu=mu,
                roughness=roughness
            )
            mdot_man_rem = mdot_man_next

    total_ch_flow = sum(mdot_channels[1:], mdot_channels[0]) * scale_factor
    Balance(
        "Total_Mass_Balance", net,
        variable=P_inlet,
        function=(total_ch_flow - mdot_half_target) / mdot_half_target
    )

    SteadyState(net).solve(solver_method="trf", ftol=1e-8, xtol=1e-8, rtol=1e-6, verbose=False)

    channel_flows = np.array([ch.value for ch in mdot_channels])
    p_static_ring = np.array([p.value for p in P_nodes])

    # Spatial static pressure variation (allows momentum regain while keeping profile smooth/uniform)
    p_static_std = float(np.std(p_static_ring))
    p_static_cov_pct = 100.0 * p_static_std / np.mean(p_static_ring)

    return {
        "mdot_distribution": channel_flows,
        "P_static_manifold": p_static_ring,
        "P_static_std": p_static_std,
        "P_static_cov_pct": p_static_cov_pct,
        "maldistribution_pct": 100.0 * (np.max(channel_flows) - np.min(channel_flows)) / np.mean(channel_flows),
        "dP_static_end_to_end": p_static_ring[0] - p_static_ring[-1],
        "P_inlet_required": P_inlet.value,
        "dP_total": P_inlet.value - P_exit
    }