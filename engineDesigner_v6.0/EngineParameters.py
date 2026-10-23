# --- 1. PROPELLANTS ---
fuel = "RP-1"
oxidizer = "LOX"

# --- 2. OPERATING CONDITIONS ---
Pc = 261.6            # psia, Chamber Pressure
Pe = 11.0              # psia, Exit Pressure (for perfect expansion)
T_target = 3597.0      # lbf, Target Actual Thrust
MR = 2.1               # Mixture Ratio
Cstar_eff = 0.85        # C* Efficiency 
Cf_eff = 0.95          # Cf Efficiency 
mdot = 7.95               #kg/s, mass flow rate

# --- 3. Coolant Conditions --- 
coolant = "RP-1"
T_cool_in = 298.15 #K, coolant inlet temperature
P_cool_out_desired = 314*6894.75729 #Desired Pressure (Pa) of Coolant before injector, 24% injector stiffness


