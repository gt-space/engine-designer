#_________________________________________________
#        TRANSONIC FLOW MODULE (TA)
#__________________________________________________

import numpy as np
import ode

class transonicFlow:
    def __init__(self, N, r_w, z_w, z_axis, gamma_arr, p_star_arr, a_star_arr, m_dot_arr, RoC, eq_solver = None):
        self.N = int(N) #num of zones
        self.r_w = float(r_w) #radial coord of wall
        self.z_w = float(z_w) #axial coord of wall
        self.z_axis = float(z_axis)# axial coord where the Mach number at r=0, z=zaxis is the same as at point rw, zw
        self.eq_solver = eq_solver

        self.gamma = np.array(gamma_arr, dtype=float) #gamma table
        self.p_star = np.array(p_star_arr, dtype=float)#P_star, provided by ODE/K or cantera
        self.a_star = np.array(a_star_arr, dtype=float)#a star, can probably get from cantera
        self.m_dot = np.array(m_dot_arr, dtype=float) #partial mass flow rates, can get from geometry or cantera
        self.RoC = float(RoC) #roc in throat
    
        self.B0 = np.zeros(self.N)
        self.B1 = np.zeros(self.N)
        self.C1 = np.zeros(self.N)
        self.C2 = np.zeros(self.N)

        self.Y = np.zeros(self.N + 1) #slipline locs, including center axis, r=0 bc axisymmetric
        self.Y[0] = 0.0
        self.Y[self.N] = self.r_w #inner and outer boundaries, wall adn 0 axis

    def trim(self, A_mat, b_vec):
        # This solves a tridiagonal matrix, we don't need all the functions they used but (5-118)
        n = len(b_vec)
        A = np.zeros((n,n))
        np.fill_diagonal(A, A_mat[:,1])        
        r = np.arange(1,n)
        A[r,r-1] = A_mat[1:,0]
        A[r-1,r] = A_mat[:-1,2]
        x = np.linalg.solve(A,b_vec)
        return x

    def bandi(self, A, B): #idt we have to do anything fancy like they did in the nasa doc, (5-105 to 5-107)
        x = np.linalg.solve(A, B)
        return x

    def sauer(self, r, x_coord, gamma): 
        #ass striated flow with constant gamma
        # Compute local flow properties and streamline angle using Sauer's method (5-116 to 5-117)
        R_prime = r + gamma/4 #what is gamma? need to find out bc this is gamma/4 
        X = np.sqrt(r/R_prime)*x_coord
        n_ind = 0
        for i in range(self.N):
            if self.Y[i]<=r<=self.Y[i+1]:
                n_ind = i
                break
        gamma_n = self.gamma[n_ind] #once again what is gamma need to find out
        
        if self.N ==1: #there's no reason for it to be really but it's in the nasa doc so I put it in
            B1sq = 2.0 / ((gamma_n + 1.0)*R_prime)
            self.B1[0] = np.sqrt(B1sq)
            self.B0[0] = -0.25/R_prime
            u_prime = 0.25*(gamma_n + 1.0)*self.B1[0]**2*r**2+ self.B0[0] + self.B1[0]*X
            v_prime = 1/16*(gamma_n + 1.0)**2*self.B1[0]**3*r**3 + 1/2*(gamma_n + 1.0)*self.B1[0]*self.B0[0]*r + 1/2*(gamma_n + 1.0)*self.B1[0]**2*r*X
        else:
            B1sq = 2.0 /((gamma_n + 1.0)*R_prime) #not sure if this is the same formula, it doesn't say but I assume so 
            self.B1[n_ind] = np.sqrt(B1sq)
            self.B0[n_ind] = -0.25/R_prime    
            c1 = self.C1[n_ind] #idk where to find these, look into
            c2 = self.C2[n_ind]
            u_prime = 0.25*(gamma_n + 1.0)*self.B1[n_ind]**2*r**2 + c1*np.log(r) + self.B0[n_ind] + self.B1[n_ind]*X
            v_prime = (1/16*(gamma_n + 1.0)**2*self.B1[n_ind]**3*r**3 + 1/2*(gamma_n + 1.0)*self.B1[n_ind]*c1*r*(np.log(r) - 1/2) + 1/2*(gamma_n + 1.0)*self.B1[n_ind]*self.B0[n_ind]*r + c2*r + (1/2*(gamma_n + 1.0)*self.B1[n_ind]**2*r + c1/r)*x_coord) #should this last one be X instead?
        u = 1.0 + u_prime
        v = v_prime
        p_over_p_star = 1.0 - gamma_n*u_prime
        theta = np.arctan(v/u)
        
        return p_over_p_star, theta, n_ind + 1
    
    def fcalc(self):
        # calculates the F functions at streamline boundaries and at the walls, (5-108 to 5-109)
        f_vals = []
        for n in range(self.N-1):
            Y_n = self.Y[n+1]
            g_n = self.gamma[n]
            g_next = self.gamma[n+1]
            B1_n = self.B1[n]
            B1_next = self.B1[n+1]
            B0_n = self.B0[n]
            B0_next = self.B0[n+1]
            C1_n = self.C1[n]
            C1_next = self.C1[n+1]
            C2_n = self.C2[n]
            C2_next = self.C2[n+1]
            P_star_n = self.p_star[n]
            P_star_next = self.p_star[n+1]

            ln_Yn = np.log(Y_n)

            F1_n = 1/2*(g_n+1.0)*B1_n**2*Y_n + C1_n/Y_n - (1/2*(g-next+1.0)*B1_next**2*Y_n + C1_next/Y_n)
            F2_n = (1/16*(g_n+1.0)**2*B1_n**3*Y_n**3 + 1/2*(g_n+1.0)*B1_n*C1_n*Y_n*(ln_Yn-1/2) + 1/2*(g_n+1.0)*B1_n*B0_n*Y_n + C2_n/Y_n) - (1/16*(g_next+1.0)**2*B1_next**3*Y_n**3 + 1/2*(g_next+1.0)*B1_next*C1_next*Y_n*(ln_Yn-1/2) + 1/2*(g_next+1.0)*B1_next*B0_next*Y_n + C2_next/Y_n)
            F3_n = g_n*B1_n*P_star_n-g_next*B1_next*P_star_next
            F4_n = (g_n*(0.25*(g_n+1.0)*B1_n**2*Y_n**2 + C1_n*ln_Yn+B0_n)-1)*P_star_n - (g_next*(0.25*(g_next+1.0)*B1_next**2*Y_n**2 + C1_next*ln_Yn + B0_next)-1)*P_star_next

            f_vals.append(F1_n)
            f_vals.append(F2_n)
            f_vals.append(F3_n)
            f_vals.append(F4_n)

        N_ind = self.N-1
        F1_wall = 1/self.RoC - (1/2*(self.gamma[N_ind]+1.0)*self.B1[N_ind]**2 + self.C1[N_ind]) #is this the radius of convergence? note to self, find out
        F2_wall = 1/16*(self.gamma[N_ind]+1.0)**2*self.B1[N_ind]**3 - 1/4*self.gamma[N_ind]
        f_vals.append(F1_wall)
        f_vals.append(F2_wall)
        return np.array(f_vals)
    
    def tran(self, b_arr, g_thermo_fn, P_c):
        # constructs initial data line for MoC by computing expansion coefs with ODE, solving slipline locs,
        # and calculating point distributions ()

        for n in range(self.N):
            # gonna need an actual temperature instead of just 2750, whatever the cantera gang says
            n_st, T_st, A_st = ode.tpcalc(self.eq_solver, 2750, P_c, b_arr, g_thermo_fn)

        A_mat = np.zeros((self.N-1, 3))
        b_vec = np.zeros(self.N-1)
        for n in range(self.N-1):
            An = self.rho_star[n] * self.a_star[n] / (self.rho_star[n+1] * self.a_star[n+1])
            kn = self.m_dot[n+1]/self.m_dot[n]
            if n>0:
                A_mat[n, 0] = -An*kn
            A_mat[n,1] = 1.0 + An*kn
            if n<self.N - 2:
               A_mat[n, 2] = -1.0
            b_vec[n] = 0.0
        self.Y[1:self.N] = self.trim(A_mat, b_vec)

        self.newt()

        initial_line_points = []

        for i in range(self.N+1):
            r_i = self.r_w * np.sin(i/self.N * np.pi/2.0)
            z_i = self.z_w + ((self.r_w-r_i) / self.r_w)**2 * (self.z_axis-self.z_w)

            p_ratio, theta, zone_index = self.sauer(r_i, z_i, self.gamma) #gamma from cantera prob
            initial_line_points.append({'r': r_i, 'z': z_i, 'P_over_P_star': p_ratio, 'theta': theta, 'zone': zone_index})
        
        return initial_line_points
        
    def get_vars(self, vars):
        # turns the one row array into the individual arrays for all B and C
        size = self.N
        self.B0 = vars[0:size]
        self.B1 = vars[size:2*size]
        self.C1 = vars[2*size:3*size]
        self.C2 = vars[3*size:4*size]

    def newt(self, eps = 1e-4, threshold = 1e-6, max_iter = 100):
        # Newton's method to solve the equations with the boundary conditions from fcalc()
        # uses the jacobian method

        vars_vec = np.concatenate([self.B0, self.B1, self.C1, self.C2])
        n_vars = len(vars_vec)
        
        for iter in range(max_iter):
            f_0 = self.fcalc()
            error_sum = np.sum(f_0**2)
            if error_sum < threshold:
                return True, iter

        num_funcs = len(f_0)
        J = np.zeros((num_funcs, n_vars))

        for j in range(n_vars):
            delta = max(abs(eps * vars_vec[j]), eps*0.001)
            vars_vec[j] += delta
            self.get_vars(vars_vec)
            f_perturbed = self.fcalc()
            J[:,j] = (f_perturbed - f_0) / delta
            vars_vec[j] -= delta
            self.get_vars(vars_vec)

        try:
            step = np.linalg.solve(J.T @ J + 1e-8 * np.eye(n_vars), -J.T @f_0)
        except np.linalg.LinAlgError:
            return False, iter if 'iter' in locals() else max_iter

        vars_vec += step
        self.get_vars(vars_vec)

        return False, max_iter

