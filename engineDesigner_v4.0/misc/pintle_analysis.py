import numpy as np

diam_pintle = 0.0233 # m

rho_LOX = 1036 # kg/m3
rho_fuel = 800 # kg/m3

P_inj_LOX = 21.24 * 100000
P_inj_fuel = 21.6 * 100000

P_c = 18 * 100000

dP_fuel = P_inj_fuel - P_c
dP_LOX = P_inj_LOX - P_c

#Primary slots
#mdot_LOX_1 = 0.9545/2 # kg/s
slot_h_1 = 0.00095 # m
slot_w_1 = 0.001016 # m
A_LOX_1 = slot_h_1 * slot_w_1
N_slots_1 = 24
A_LOX_tot_1 = A_LOX_1 * N_slots_1
Cd_LOX_1 = 0.7638

# Secondary slots
#mdot_LOX_2 = 0.9545/2 # kg/s
slot_h_2 = 0.00095 # m
slot_w_2 = 0.00102 # m
A_LOX_2 = slot_h_1 * slot_w_1
N_slots_2 = 1
A_LOX_tot_2 = A_LOX_1 * N_slots_2
Cd_LOX_2 = 0.7638

# Fuel

# mdot_fuel = 0.5303 #kg/s
Cd_fuel = 0.76
A_fuel = 0.00002907 # m^2
ann_gap = 0.0003901 # m


mdot_fuel = Cd_fuel * A_fuel * np.sqrt(rho_fuel * 2 * dP_fuel)
vel_fuel = Cd_fuel * np.sqrt(2*dP_fuel/rho_fuel)

mdot_LOX_1 = Cd_LOX_1 * A_LOX_tot_1 * np.sqrt(2*rho_LOX*dP_LOX)
vel_LOX_1 = Cd_LOX_1 * np.sqrt(2*dP_LOX/rho_LOX)
mdot_LOX_2 = Cd_LOX_2 * A_LOX_tot_2 * np.sqrt(2*rho_LOX*dP_LOX)
vel_LOX_2 = Cd_LOX_2 * np.sqrt(2*dP_LOX/rho_LOX)

print(vel_fuel) # m/s
print(vel_LOX_1) # m/s
print(vel_LOX_2) # m/s