#____________________________________________________________
#        ONE-DIMENSIONAL EQUILIBRIUM (ODE)
# ____________________________________________________________

import numpy as np


class CPHS: #just getting the CP, H, S values and entropy using curve fits, we can probably use cantera for this 
    def evaluate(a, T):
        # a = array of 
        Cp_R = a[0]+ a[1]*T + a[2]*T**2 +a[3]*T**3 + a[4]*T**4
        H_RT = a[0] + (a[1]*T)/2.0 + (a[2]*T**2)/3.0+ (a[3]*T**3)/4.0 + (a[4]*T**4)/5.0 + a[5]/T
        S_R  = a[0]*np.log(T)+ a[1]*T + (a[2]*T**2)/2.0 + (a[3]*T**3)/3.0 + (a[4]*T**4)/4.0 + a[6]
        G_RT = H_RT - S_R
        return Cp_R, H_RT, S_R, G_RT

# Solve w equilibrium (5.2.5, 5-31 in NASA SN 63) and not
class SOLVE:
    def __init__(self, species, gas_indices, cond_indices, stoich, P0=1.0, R=8.31446):
        self.species_name = species # str names of species
        self.gas_idx = np.array(gas_indices, dtype=int)# ind of gaseous species
        self.cond_idx = np.array(cond_indices, dtype=int)# indices of condensed species
        self.m = len(self.gas_idx)
        self.n = len(species)
        self.l = len(stoich)
        self.a = np.array(stoich, dtype=float)     # [numElements, numSpecies] of a_{k,i}, stoich coefficients
        self.P0 = P0 # stag pressure
        self.R = R # universal gas constant in J/mol-K
    def lambda_calc(self, d_ln_T, d_ln_A, d_ln_n, n_curr): # empirical convergence parameter controls correction size in N-R method, 0<lambda<=1
        lambda1_poss = [abs(d_ln_T)/2, abs(d_ln_A)/2]
        for j in range(self.m):
            if d_ln_n[j] > 0:
                ln_ratio = np.log(n_curr[self.gas_idx[j]] / self.P0)
                if ln_ratio > -18.5:
                    lambda1_poss.append(d_ln_n[j])
        max_value = max(lambda1_poss)
        if max_value > 0:
            lambda_1 = 1.0 / max_value
        else:
            lambda_1 = 1.0
        lambda2_poss = []
        for j in range(self.m):
            if d_ln_n[j] > 0:
                ln_ratio = np.log(n_curr[self.gas_idx[j]]/self.P0)
                if ln_ratio <= -18.5:
                    lam2 = (np.log(self.P0) - 9.212 - np.log(n_curr[self.gas_idx[j]]))/ d_ln_n[j]
                    lambda2_poss.append(lam2)
        if lambda2_poss:
            lambda_2 = min(lambda2_poss)
        else:
            lambda2 = 1.0
        
        return min(lambda_1, lambda_2, 1.0)
    
    def converged(self, d_ln_n_gas, d_n_cond, n_tot):          # makes sure the convergence critera on 5-38 are met
        sum_n = np.sum(n_tot)
        for j_idx, gas_i in enumerate(self.cond_idx):
            check_val = (n_tot[gas_i] / sum_n) * abs(d_ln_n_gas[j_idx])
            if check_val >= 0.5e-5:
                return False
        for j_idx, cond_i in enumerate(self.cond_idx):
            check_val = abs(d_n_cond[j_idx]) / sum_n
            if check_val >= 0.5e-5:
                return False
        return True
    
    def iter_matrix_sys(self, n, T, A, b_elements, g_thermo, mode='TP'):
        # n = mole vector for all species
        # T = temperature
        # A = free-energy parameter
        # b_elements = vector of element abundances
        # g_thermo = dictionary with H_i, S_i, C_i for all species
        # mode =string with inputs type (TP, HP, SP, need two variables, should be getting TP from cantera gang)
        # 0 to l-1 are element variables delta ln u_k
        # l to l + num_cond -1 are moles of condensed species Delta n_j
        # l+num_cond is scale variable -Delta lnA
        # l+num_cond+1 is temperature variable Delta lnT
        l = self.l
        num_cond = len(self.cond_idx)
        mat_dim = l + num_cond + 2
        M = np.zeros((mat_dim, mat_dim))
        right_side = np.zeros(mat_dim)
        
        r = np.zeros((l,l))
        for k in range(l):
            for j in range(l):
                r[k,j] = np.sum(self.a[k, self.gas_idx] * self.a[j, self.gas_idx] * n[self.gas_idx])
        
        for k in range(l):
            M[k, 0:l] = r[k, :]
            
            for j_idx, cond_i in enumerate(self.cond_idx):
                M[k, l+j_idx] = self.a[k, cond_i]
            
            M[k,l+num_cond] = -np.sum(self.a[k, self.gas_idx] * n[self.gas_idx])

            H_gas = g_thermo['H_i'][self.gas_idx]
            M[k, l+num_cond+1] = np.sum(self.a[k, self.gas_idx] * n[self.gas_idx] * H_gas)

            a_all = np.sum(self.a[k,:] * n)
            right_side[k] = A * b_elements[k] + a_all

        row_t = l + num_cond + 1
        if mode == 'TP':
            M[row_t, row_t] = 1.0
            right_side[row_t] = 0.0
        elif mode =='HP':
            H_all = g_thermo['H_i']
            M[row_t, 0:l] = [np.sum(self.a[k, self.gas_idx] * H_gas * n[self.gas_idx]) for k in range(l)]
            for j_idx, cond_i in enumerate(self.cond_idx):
                M[row_t, l + j_idx] = H_all[cond_i]
            M[row_t, l + num_cond] = -np.sum(H_gas * n[self.gas_idx])
            C_gas = g_thermo['C_i'][self.gas_idx]
            M[row_t, l+num_cond+1] = np.sum(C_gas * n[self.gas_idx]) + np.sum(H_gas**2 * n[self.gas_idx])
            right_side[row_t] = g_thermo['H_assigned'] - np.sum(H_all * n)
        return M, right_side
    
    def newton_raphson_step(self, n, T, A, b_elements, g_thermo, mode='TP'):
        # make the system Mx = right_side
        M, right_side = self.iter_matrix_sys(n, T, A, b_elements, g_thermo, mode)
        corr_vec = np.linalg.solve(M, right_side)
        d_ln_u = corr_vec[0:self.l]
        d_n_cond = corr_vec[self.l:self.l + len(self.cond_idx)]
        d_ln_A = -corr_vec[self.l + len(self.cond_idx)]
        d_ln_T = corr_vec[self.l + len(self.cond_idx) + 1]
        # find Delta ln n_i for all gases from 5-31
        d_ln_n_gas = np.zeros(self.m)

        for i_idx, gas_i in enumerate(self.gas_idx):
            S_i = g_thermo['S_i'][gas_i]
            H_i = g_thermo['H_i'][gas_i]
            sum_u = np.sum(self.a[:, gas_i] * d_ln_u)
            d_ln_n_gas[i_idx] = -S_i + sum_u + H_i * d_ln_T
        #calculate scaling parameter (5-37)
        lam = self.lambda_calc(d_ln_T, d_ln_A, d_ln_n_gas, n)

        # apply corrections
        n_new = np.copy(n)
        for i_idx, gas_i in enumerate(self.gas_idx):
            n_new[gas_i] = np.exp(np.log(n[gas_i]) + lam * d_ln_n_gas[i_idx])
        
        for j_idx, cond_i in enumerate(self.cond_idx):
            n_new[cond_i] = n[cond_i] + lam* d_n_cond[j_idx]
        
        A_new = np.exp(np.log(A) + lam*d_ln_A)
        T_new = np.exp(np.log(T) + lam* d_ln_T)

        conv = self.converged(d_ln_n_gas, d_n_cond, n_new)

        return n_new, T_new, A_new, conv
    
    def thermo_derivatives(self, n, T, A, g_thermo):
        M, _ = self.iter_matrix_sys(n, T, A, np.zeros(self.l), g_thermo, mode = 'TP')
        right_side_53 = np.zeros(M.shape[0])
        for k in range(self.l):
            right_side_53[k] = -np.sum(g_thermo['H_i'][self.gas_idx] * self.a[k, self.gas_idx] * n[self.gas_idx])

        sol_53 = np.linalg.solve(M, right_side_53)
        d_ln_u_dT = sol_53[0:self.l]
        d_n_cond_dT = sol_53[self.l:self.l + len(self.cond_idx)]
        d_ln_A_dT = sol_53[self.l + len(self.cond_idx)]

        sum1 = 0.0
        for k in range(self.l):
            for i_idx, gas_i in enumerate(self.gas_idx):
                sum1 += self.a[k, gas_i] * n[gas_i] * d_ln_u_dT[k]
        
        sum2 = np.sum(g_thermo['H_i'][self.cond_idx] * d_n_cond_dT)
        sum3 = np.sum(g_thermo['H_i'][self.gas_idx] * n[self.gas_idx] * -d_ln_A_dT)
        sum4 = np.sum(g_thermo['C_i'][self.gas_idx] * n[self.gas_idx])
        sum5 = np.sum(g_thermo['H_i'][self.gas_idx] * n[self.gas_idx])

        Cp_eq = self.R * (sum1 + sum2 + sum3+sum4 + sum5)
        d_ln_Mw_d_ln_T = d_ln_A_dT
        
        return Cp_eq, d_ln_Mw_d_ln_T
    def frozen(n_frozen, S_comb, P_exit, T_initial, g_thermo_fn, R = 8.314463, max_iter=5000, tol=1e-5):
        # find temperature and thermodynamic props for frozen composition expanded to exit pressure, S = combustion entropy, returns dictionary
        # n_frozen     = fixed comp mole array frozen @ combustion conditions
        # S_comb       = entropy at combustion conditions
        # P_exit       = exit pressure
        # T_initial    = initial estimate for T
        # g_thermo_fn  = dictionary with S, C, Mw, m_gas
        # R            = universal gas constant
        # max_iter     = when to stop if not stopping
        # tol          = convergence tolerance
        T = float(T_initial)
        m = len(n_frozen)

        for iter in range(max_iter):
            thermo = g_thermo_fn(T, P_exit)
            S_i = thermo['S_i']
            C_i = thermo['C_i']

            S_tot = np.sum(S_i * n_frozen) #total entropy
            
            C_sum = np.sum(C_i * n_frozen)
            d_ln_T = (S_comb - S_tot) / C_sum

            T_new = np.exp(np.log(T) + d_ln_T)

            if abs(T_new-T) < tol:
                T = T_new
                break
            T = T_new
        
        thermo = g_thermo_fn(T, P_exit)
        C_i = thermo['C_i']
        Mw_i = thermo['Mw_i']

        m_gas = thermo.get('m_gas', len(n_frozen))

        x_i = n_frozen
        sum_x_gas = np.sum(x_i[:m_gas])
        Mw_c = np.sum(x_i[:m_gas] * Mw_i[:m_gas]) / sum_x_gas if sum_x_gas >0 else np.sum(x_i * Mw_i)

        Cp = R * np.sum(C_i*x_i)/Mw_c*sum_x_gas

        d_ln_Mw_d_ln_T = 0.0
        d_ln_Mw_d_ln_P = 0.0

        gamma = (Cp * Mw_c)/(Cp*Mw_c-R)

        return {'T': T, 'Cp': Cp, 'd_ln_Mw_d_ln_T': d_ln_Mw_d_ln_T, 'd_ln_Mw_d_ln_P': d_ln_Mw_d_ln_P, 'gamma': gamma, 'Mw_c': Mw_c}
    
    def frozen_throat(T_c, P_c, gamma_c, h_c, Mw_c, g_thermo_fn, R = 8.314463, max_iter=50, tol=1e-5):
        # gives throat pressure and temperature, frozen
        # T_c = chamber temperature
        # P_c = chamber pressure
        # gamma_c = gamma at combustion conditions
        # h_c = specific enthalpy, combustion conditions
        # Mw_c = mixture molecular weight, combustion conditions
        # g_thermo_fn = same as before

        P_ratio_est = ((gamma_c + 1.0)/2.0) ** (gamma_c/(gamma_c-1.0))
        T_t = T_c * (2.0/(1.0+gamma_c))
        P_t = P_c / P_ratio_est

        for k in range(max_iter):
            thermo_t = g_thermo_fn(T_t, P_t)
            h_t = thermo_t['h_total']
            h_diff = abs(h_c-h_t)/abs(h_c)
            if h_diff <= tol:
                break
            P_ratio_est = P_ratio_est / (1.0 + (2.0 *Mw_c*(h_c-h_t))/((gamma_c+1.0)*R*T_t))
            P_t = P_c / P_ratio_est
            T_t = T_c * (2.0/(1.0+gamma_c))
        
        return P_t, T_t

class SpeciesSelect:
    def __init__(self, prop_db):
        # I haven't looked at Tyson's properties db yet but this may connect there
        self.db = prop_db

    # figure this out

def tpcalc(eq_solver, T, P, b_arr, g_thermo, max_iter = 5000):
    # eq_solver = instance of  Solve
    # T = temp (K)
    # P = pressure #don't think it matters what unit
    # b_arr = relative abundances
    # g_thermo = dictionary with H_i, S_i, C_i for all species
    eq_solver.P0 = P
    n = np.ones(eq_solver.n) * 0.1/eq_solver.m
    A = 1.0
    for _ in range(max_iter):
        g_th = g_thermo(T)
        n, T, A, conv = eq_solver.newton_raphson_step(n, T, A, b_arr, g_th, mode='TP')
        if conv:
            break

    return n, T, A

def rocket_out(h_c, h, P_c, A_w_t, P, A_w, gamma, T, Mw):
    #I think these are the right units but we I need to check because of the scaling constants that the nasa people used 
    # h_c = combustion chamber h (kcal/kg)
    # h = exit h (cal/g)
    # P_c = comb chamber pressure (psi)
    # A_w_t = throat area per mass flow rate (s*in^2/lbm)
    # P = exit pressure (psi)
    # A_w = exit area per mass flow rate (s*in^2/lbm)
    # gamma = specific heat ratio
    # T = local temp (K)
    # Mw = mol weight (kg/kmol)
    g_c = 32.174 # ft/s^2
    I = 294.94 * np.sqrt((h_c-h)/1000.0)
    I_vac = I + P * A_w
    c_star = g_c * P_c * A_w_t
    C_F = g_c*I/c_star
    M = I/np.sqrt(86.4579*gamma*T/Mw)

    return {'I': I, 'I_vac': I_vac, 'c_star': c_star, 'C_F': C_F, 'M': M}

def rocket(eq_solver, P_c, P_e, b_arr, g_thermo, mode='eq')
    n_c, T_c, _ = tpcalc(eq_solver, 3000.0, P_c, b_arr, g_thermo)
    g_th_c = g_thermo(T_c)
    h_c = np.sum(n_c*g_th_c['H_i'])
    S_c = np.sum(n_c*g_th_c['S_i'])
    if mode == 'frozen':
        f = eq_solver.frozen(n_frozen = n_c, S_c=S_c, P_exit = P_e, T_initial = T_c*0.8, g_thermo = g_thermo)
        T_e = f['T']
        gamma_e = f['gamma']
        Mw_e = f['Mw_c']
        g_th_e = g_thermo(T_e, P_e)
        h_e = np.sum(n_c*g_th_e['H_i'])
    elif mode=='eq':
        n_e, T_e, A_e = tpcalc(eq_solver, 1500.0,P_e, b_arr, g_thermo)
        g_th_e = g_thermo(T_e, P_e)
        h_e = np.sum(n_e* g_th_e['H_i'])
        sum_x_gas = np.sum(n_e[:eq_solver.m])
        Mw_i = g_th_e['Mw_i']
        Mw_e = []
        np.sum(n_e[:eq_solver.m]*Mw_i[:eq_solver.m])/sum_x_gas if sum_x_gas>0 else np.sum(n_e*Mw_i)
        Cp_eq, _ = eq_solver.thermo_derivatives(n_e, T_e, A_e, g_th_e)
        gamma_e = (Cp_eq * Mw_e)/(Cp_eq*Mw_e-eq_solver.R)
    # add something to break if the mode isn't right
    perf_params = rocket_out(h_c = h_c, h = h_e, P_c = P_c, A_w_t = [], P = P_e, A_w = [], gamma=gamma_e, T=T_e, Mw = Mw_e)
    perf_params['T_e'] = T_e
    perf_params['gamma_e'] = gamma_e
    perf_params['Mw_e'] = Mw_e
    
    return perf_params