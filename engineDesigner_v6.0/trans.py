#_________________________________________________
#        TRANSONIC FLOW MODULE (TA)
#__________________________________________________

import numpy as np
import ode

class TransonicFlow:
    def __init__(self, N, r_w, z_w, z_axis, gamma_arr, p_star_arr, a_star_arr, m_dot_arr, RoC):
        self.N = int(N) #num of zones
        self.r_w = float(r_w) #radial coord of wall
        self.z_w = float(z_w) #axial coord of wall
        self.z_axis = float(z_axis)# axial coord of axis
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
        n = len(b_vec)
        A = np.zeros((n,n))
        np.fill_diagonal(A, A_mat[:,1])        
        r = np.arange(1,n)
        A[r,r-1] = A_mat[1:,0]
        A[r-1,r] = A_mat[:-1,2]
        x = np.linalg.solve(A,b_vec)
        return x

    def bandi(self, A, B): #idt we have to do anything fancy like they did in the nasa doc
        x = np.linalg.solve(A, B)
        return x

    def sauer(self, r, x_coord, gamma): #ass striated flow with constant gamma
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