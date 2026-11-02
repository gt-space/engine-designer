# --- 1. Propellents ---
fuel = "RP-1"
oxidizer = "LOX"

# --- 2. Operating Conditions ---
Pc = 261.6            # psia, Chamber Pressure
Pe = 11.0              # psia, Exit Pressure (for perfect expansion)
T_target = 3597.0      # lbf, Target Actual Thrust
MR = 2.1               # Mixture Ratio
Cstar_eff = 0.85        # C* Efficiency 
Cf_eff = 0.95          # Cf Efficiency 
mdot = 7.95               #kg/s, mass flow rate
D_c = 6.5               # in, Diameter of the Injector
Rtd_Rt = 0.382          # Fixed Rao Optimal Constant

# --- 3. Coolant Conditions --- 
coolant = "RP-1"
T_cool_in = 298.15 #K, coolant inlet temperature
P_cool_out_desired = 314*6894.75729 #Desired Pressure (Pa) of Coolant before injector, 24% injector stiffness

# --- 4. Other Constants ---
wall_material = "42"
channel_roughness = 10.0e-6
tw_default = 1.0e-3