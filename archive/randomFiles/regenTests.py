# Look into this. Doesn't seem to match the H&H formula
h_c = 0.021 * (Re_c**0.8) * (Pr_c**0.4) * (0.64 + 0.36 * (T_cb/T_wc)) * (cond_c/D_hyd) #coolant side
m = math.sqrt((h_c * perim_conf)/(cond_w * A_conf)) #fin efficiency parameter
eta_fin = math.tanh(m * L_cor)/(m * L_cor) #fin efficiency
eta_tot = 1 - (((numChannels * A_fin)/A_totc) * (1 - eta_fin)) #overall coolant side heat transfer efficiency
q_ce = h_c * A_totc * eta_tot * (T_wc - T_cb) #Coolant side heat transfer estimate
# print(h_c)
if abs((q_ce - q_ge)/q_ge) < 0.001:
    break
elif q_ce > q_ge:
    T_wc -= (q_ce - q_ge)/q_ce
else:
    T_wc += (q_ge - q_ce)/q_ge
