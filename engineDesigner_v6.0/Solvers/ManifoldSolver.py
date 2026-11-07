import numpy as np
from fullflow import Network, State, Component, Balance, SteadyState, DarcyWeisbach
from Functions import PressureDrop as pd

#Single feed pipe splitting into two symmetric 180-deg manifold streams at split angle alpha
class ManifoldInletSplit(Component):

    def __init__(self, name, network, mdot_pipe, P_pipe, P_man_inlet, 
                 rho, A_pipe, A_man_inlet, split_angle_deg=90.0):
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
        # Scale residual by 1e5 Pa for numerical conditioning
        self.residual = ((self.P_pipe.value - dP_static) - self.P_man_inlet.value) / 1e5

    @property
    def balances(self):
        return [(self.P_man_inlet, self.residual)]

#Circumferential manifold segment between channel i and channel i+1.
# Balances wall friction against deceleration momentum static pressure regain.
class ManifoldRingSegment(Component):

    def __init__(self, name, network, P_up, P_down, mdot_up, mdot_down, 
                 A_up, A_down, Dh_up, Dh_down, dL, rho, mu, roughness):
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
        # P_down = P_up + dP_static_change
        self.p_seg_residual = ((self.P_up.value + dP_static_change) - self.P_down.value) / 1e5

    @property
    def balances(self):
        return [(self.P_down, self.p_seg_residual)]

#90-degree branch takeoff + area contraction from manifold station i into cooling channel i.
class ManifoldBranchNode(Component):

    def __init__(self, name, network, P_man_local, P_ch_entry, mdot_man_in, 
                 mdot_ch, A_man, A_ch, rho):
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
    
    N_half = channel_geom['N_channels'] // 2
    mdot_total = engine_params['mdot_total']
    mdot_half_target = 0.5 * mdot_total
    rho = engine_params['rho']
    mu = engine_params['mu']
    roughness = engine_params['roughness']
    P_exit = engine_params['P_collector_exit']
    
    # Arc length between adjacent cooling channels along the manifold ring
    d_theta = np.pi / N_half
    dL_segment = channel_geom['R_manifold'] * d_theta

    # Unknown inlet pressure (varied by FullFlow Balance so sum(mdot_ch) == mdot_half_target)
    P_inlet = State(engine_params['P_inlet_guess'], bounds=(P_exit + 1e4, 5e8), keep_feasible=True)
    P_man_0 = State(engine_params['P_inlet_guess'] * 0.98, bounds=(P_exit + 1e3, 5e8), keep_feasible=True)
    
    # 1. Pipe-to-Manifold Wye/Tee Split
    ManifoldInletSplit(
        "Inlet_Split", net,
        mdot_pipe=mdot_total,
        P_pipe=P_inlet,
        P_man_inlet=P_man_0,
        rho=rho,
        A_pipe=engine_params['A_feed_pipe'],
        A_man_inlet=A_man_profile[0],
        split_angle_deg=inlet_angle_deg
    )
    
    P_nodes = [P_man_0]
    mdot_channels = []
    mdot_man_rem = State(mdot_half_target)  # Remaining circumferential flow in manifold
    
    # 2. March around the 180-deg half-ring
    for i in range(N_half):
        P_ch_in = State(P_man_0.value * 0.95, bounds=(P_exit + 100.0, 5e8), keep_feasible=True)
        mdot_ch = State(mdot_half_target / N_half, bounds=(1e-6, mdot_half_target), keep_feasible=True)
        mdot_channels.append(mdot_ch)
        
        # 90-deg turn + contraction into cooling channel i
        ManifoldBranchNode(
            f"Branch_{i}", net,
            P_man_local=P_nodes[i],
            P_ch_entry=P_ch_in,
            mdot_man_in=mdot_man_rem,
            mdot_ch=mdot_ch,
            A_man=A_man_profile[i],
            A_ch=channel_geom['A_ch_inlet'],
            rho=rho
        )
        
        # Downstream cooling channel resistance (solves for mdot_ch)
        DarcyWeisbach(
            f"CoolingChannel_{i}", net,
            mass_flow=mdot_ch,
            upstream_pressure=P_ch_in,
            downstream_pressure=P_exit,
            length=channel_geom['L_channel'],
            hydraulic_diameter=channel_geom['Dh_channel'],
            density=rho,
            cross_sectional_area=channel_geom['A_channel_avg'],
            friction_factor=channel_geom.get('f_channel', 0.02)
        )
        
        # Circumferential segment from station i to i+1
        if i < N_half - 1:
            mdot_man_next = State.maximum(mdot_man_rem - mdot_ch, 1e-6)
            P_next = State(P_nodes[i].value, bounds=(P_exit + 1e3, 5e8), keep_feasible=True)
            P_nodes.append(P_next)
            
            ManifoldRingSegment(
                f"RingSeg_{i}", net,
                P_up=P_nodes[i],
                P_down=P_next,
                mdot_up=mdot_man_rem,
                mdot_down=mdot_man_next,
                A_up=A_man_profile[i],
                A_down=A_man_profile[i+1],
                Dh_up=Dh_man_profile[i],
                Dh_down=Dh_man_profile[i+1],
                dL=dL_segment,
                rho=rho,
                mu=mu,
                roughness=roughness
            )
            mdot_man_rem = mdot_man_next

    # 3. Global Mass Conservation Closure: sum of channels must equal half of total flow
    total_ch_flow = sum(mdot_channels[1:], mdot_channels[0])
    Balance(
        "Total_Mass_Balance", net,
        variable=P_inlet,
        function=(total_ch_flow - mdot_half_target) / mdot_half_target
    )
    
    # Solve steady-state network
    SteadyState(net).solve(solver_method="trf", ftol=1e-7, xtol=1e-7, rtol=1e-2, verbose=False)
    
    channel_flows = np.array([ch.value for ch in mdot_channels])
    p_static_ring = np.array([p.value for p in P_nodes])
    
    return {
        "mdot_distribution": channel_flows,
        "P_static_manifold": p_static_ring,
        "maldistribution_pct": 100.0 * (np.max(channel_flows) - np.min(channel_flows)) / np.mean(channel_flows),
        "dP_static_end_to_end": p_static_ring[0] - p_static_ring[-1],
        "P_inlet_required": P_inlet.value,
        "dP_total": P_inlet.value - P_exit
    }