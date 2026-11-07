from dataclasses import dataclass, asdict
import json
import numpy as np
import matprotlib as mp


def load_matprotlib_material(mat_input):
    """Helper to resolve shorthand codes or full names into a matprotlib Material object."""
    if hasattr(mat_input, "get"):
        return mat_input
    key = str(mat_input).strip().lower()
    shortcuts = {
        "42": lambda: mp.GrCop("42"),
        "grcop-42": lambda: mp.GrCop("42"),
        "84": lambda: mp.GrCop("84"),
        "grcop-84": lambda: mp.GrCop("84"),
        "718": lambda: mp.Inconel("718"),
        "inconel 718": lambda: mp.Inconel("718"),
        "625": lambda: mp.Inconel("625"),
        "inconel 625": lambda: mp.Inconel("625"),
        "c101": lambda: mp.Copper("C101"),
        "6061": lambda: mp.Aluminum("6061"),
        "304": lambda: mp.Stainless("304"),
        "316": lambda: mp.Stainless("316"),
    }
    if key in shortcuts:
        return shortcuts[key]()
    return mp.get_material(str(mat_input))


@dataclass
class EngineConfig:
    # --- 1. Propellants ---
    fuel: str = "RP-1"
    oxidizer: str = "LOX"

    # --- 2. Operating Conditions ---
    Pc: float = 261.6              # psia, Chamber Pressure
    Pe: float = 11.0               # psia, Exit Pressure
    T_target: float = 3597.0       # lbf, Target Actual Thrust
    MR: float = 2.1                # Mixture Ratio (O/F)
    Cstar_eff: float = 0.85        # C* Efficiency
    Cf_eff: float = 0.95           # Cf Efficiency
    mdot: float = 7.95             # kg/s, Total Mass Flow Rate
    D_c: float = 6.5               # in, Chamber / Injector Diameter
    Rtd_Rt: float = 0.382          # Fixed Rao Downstream Throat Fillet Ratio

    # --- 3. Coolant Conditions ---
    coolant: str = "RP-1"
    T_cool_in: float = 298.15                          # K, Coolant Inlet Temperature
    P_cool_out_desired: float = 314.0 * 6894.75729     # Pa, Desired Channel Exit Pressure
    d_feed_pipe: float = 0.75 * 0.0254                 # m, Single Inlet Feed Pipe ID

    # --- 4. Liner & Outer Shell Material / Structural Limits ---
    wall_material: str = "42"              # Legacy alias for GrCop-42
    liner_material_name: str = "GrCop-42"  # Inner hot-gas liner material in matprotlib
    shell_material_name: str = "GrCop-42" # Outer closeout shell material in matprotlib
    channel_roughness: float = 10.0e-6     # m, Cooling channel surface roughness
    tw_default: float = 1.0e-3             # m, Baseline initial liner thickness (1.0 mm)
    t_liner_min_mfg: float = 0.75e-3       # m, Min manufacturable liner thickness (0.75 mm)
    t_shell_min_mfg: float = 1.50e-3       # m, Min manufacturable outer shell thickness (1.50 mm)
    SF_structural: float = 2.0             # Yield safety factor for pressure/structural loads
    SF_thermal: float = 1.5                # Safety factor for thermal ratcheting / expansion

    @property
    def Pc_Pa(self) -> float:
        return self.Pc * 6894.75729

    @property
    def mdot_fuel(self) -> float:
        return self.mdot / (1.0 + self.MR)

    @property
    def A_feed_pipe(self) -> float:
        return 0.25 * np.pi * self.d_feed_pipe**2

    @property
    def liner_material(self):
        return load_matprotlib_material(self.liner_material_name)

    @property
    def shell_material(self):
        return load_matprotlib_material(self.shell_material_name)

    def save_json(self, filepath: str):
        with open(filepath, "w") as f:
            json.dump(asdict(self), f, indent=4)


# Instantiate active configuration
ENGINE_CONFIG = EngineConfig()

# =========================================================================
# Backward-Compatible Module-Level Exports 
# =========================================================================
fuel = ENGINE_CONFIG.fuel
oxidizer = ENGINE_CONFIG.oxidizer
Pc = ENGINE_CONFIG.Pc
Pe = ENGINE_CONFIG.Pe
T_target = ENGINE_CONFIG.T_target
MR = ENGINE_CONFIG.MR
Cstar_eff = ENGINE_CONFIG.Cstar_eff
Cf_eff = ENGINE_CONFIG.Cf_eff
mdot = ENGINE_CONFIG.mdot
D_c = ENGINE_CONFIG.D_c
Rtd_Rt = ENGINE_CONFIG.Rtd_Rt

coolant = ENGINE_CONFIG.coolant
T_cool_in = ENGINE_CONFIG.T_cool_in
P_cool_out_desired = ENGINE_CONFIG.P_cool_out_desired
d_feed_pipe = ENGINE_CONFIG.d_feed_pipe

wall_material = ENGINE_CONFIG.wall_material
liner_material_name = ENGINE_CONFIG.liner_material_name
shell_material_name = ENGINE_CONFIG.shell_material_name
channel_roughness = ENGINE_CONFIG.channel_roughness
tw_default = ENGINE_CONFIG.tw_default
t_liner_min_mfg = ENGINE_CONFIG.t_liner_min_mfg
t_shell_min_mfg = ENGINE_CONFIG.t_shell_min_mfg
SF_structural = ENGINE_CONFIG.SF_structural
SF_thermal = ENGINE_CONFIG.SF_thermal