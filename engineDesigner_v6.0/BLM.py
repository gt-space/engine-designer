"""
This is the current state of the Boundary Layer Module (BLM) as part of the TDK code

Basic Logic:
Get flow and geometry input -> Divide nozzle into N axial grid points and J points in y-direction
-> Solve boundary layer equations at each x point, march downstream -> Can compute BL thickness, heat transfer, ... from solved BL equations

Implemented:
    - Basic BL equation solver 
    - Downstram Marching
    - Esentially the code can already be run with dummy physics and calculate BL velocities and thickness

Still Missing:
    - Most of the real physics (Is done by computing coefficients)
    - Turbulence model
    - Physical grid, current grid has no physicall representation

Felix Lindner
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass


@dataclass
class Grid:
    xi: np.array
    k: np.array
    eta: np.array
    h: np.array

@dataclass
class BL:
    f: np.array
    u: np.array
    v: np.array
    g: np.array
    p: np.array

@dataclass
class P_coeff:
    m1: np.array
    m2: np.array
    m3: np.array

g_w=1
def get_input():
    max_iter = 10   #Max number of iterations per x location
    return max_iter

def build_grid(xi_max, NX, eta_max, NY):
    #Create Grid in xi and eta direction
    xi = np.linspace(0, xi_max, NX)
    k = np.diff(xi)
    eta = np.linspace(0, eta_max, NY)
    h = np.diff(eta)

    grid = Grid(xi=xi, k=k, eta=eta, h=h)
    return grid

def calc_pressure_coefficients(grid):
    xi = grid.xi
    N = len(xi)

    #Allocate
    m1 = np.zeros(N)
    m2 = np.zeros(N)
    m3 = np.zeros(N)    #m3 only for suction or blowing, =0 if neither.
    m4 = np.zeros(N)

    #Compute gradients using central differences
    due_dxi = np.zeros(N)
    drhoemue_dxi = np.zeros(N)

    for i in range(1, N-1):
        due_dxi[i] = (ue[i+1] - ue[i-1]) / (xi[i+1] - xi[i-1])
        drhoemue_dxi[i] = (rhoe[i+1] * mue[i+1] - rhoe[i-1] * mue[i-1]) / (xi[i+1] - xi[i-1])

    # Forward/backward difference at boundaries
    due_dxi[0] = (ue[1] - ue[0]) / (xi[1] - xi[0])
    drhoemue_dxi[0] = (rhoe[1] * mue[1] - rhoe[0] * mue[0]) / (xi[1] - xi[0])
    due_dxi[-1] = (ue[-1] - ue[-2]) / (xi[-1] - xi[-2])
    drhoemue_dxi[-1] = (rhoe[-1] * mue[-1] - rhoe[-2] * mue[-2]) / (xi[-1] - xi[-2])

    #Compute coefficients
    for i in range(N):
        m2[i] = (xi[i] / ue[i]) * due_dxi[i]
        m4[i] = (xi[i] / (mue[i] * rhoe[i])) * drhoemue_dxi[i]
        m1[i] = (1 + m2[i] + m4[i])/2

    return P_coeff(m1=m1, m2=m2, m3=m3)
    




def compute_local_coefficients(eta, f, u, v, g, p, xi):
    J = len(eta)

    # Allocate
    b  = np.zeros(J)
    c  = np.zeros(J)
    d  = np.zeros(J)
    e  = np.zeros(J)

    for j in range(J):

        # Density ratio
        c[j] = 1.0 / (1.0 + 0.2 * (g[j] - 1))

        # Viscosity / energy coefficients
        b[j] = (1 + 0.1 * eta[j])
        d = b.copy()
        e = b.copy()

    return b, c, d, e


def initial_bl(grid, max_iter, pressure):
    eta = grid.eta
    h = grid.h
    J = len(eta)

    #Initial Guess
    f = eta.copy()
    u = eta / np.max(eta)
    v = np.zeros(len(eta))
    g = np.zeros(len(eta))
    p = np.zeros(len(eta))

    #Pressure Coefficients
    m1 = pressure.m1[0]
    m2 = pressure.m2[0]
    m3 = pressure.m3[0]

    for it in range(1, max_iter):
        #Create matrix A and vector b for system of equations
        size  = 5 * J
        A_mat = np.zeros((size, size))
        b_vec = np.zeros(size)

        #First two rows for boundary conditions f(0) = 0 and u(0) = 0
        A_mat[0, 0] = 1.0
        b_vec[0]    = 0.0 - f[0]

        A_mat[1, 1] = 1.0
        b_vec[1]    = 0.0 - u[0]

        #g(0) = gw
        A_mat[2, 3] = 1.0
        b_vec[2] = g_w - g[0]

        b, c, d, e = compute_local_coefficients(eta, f, u, v, g, p, xi=0)


        for j in range(0, J-1):
            hj = h[j]

            #Eta Derivatives
            df_deta = (f[j+1] - f[j])/hj
            du_deta = (u[j+1] - u[j])/hj
            dg_deta = (g[j+1] - g[j])/hj
            dbv_deta = (b[j+1]*v[j+1] - b[j]*v[j])/hj
            dep_deta = (e[j+1]*p[j+1] - e[j]*p[j])/hj
            dduv_deta = (d[j+1]*u[j+1]*v[j+1] - d[j]*u[j]*v[j])/hj

            #Mid Points (j+1)
            f_mid = (f[j] + f[j+1])/2
            u_mid = (u[j] + u[j+1])/2
            v_mid = (v[j] + v[j+1])/2
            b_mid = (b[j] + b[j+1])/2
            c_mid = (c[j] + c[j+1])/2
            p_mid = (p[j] + p[j+1])/2
            g_mid = (g[j] + g[j+1])/2

            #Residuals
            R1 = df_deta - u_mid
            R2 = du_deta - v_mid
            R3 = dbv_deta + m1 * f_mid * v_mid + m2 * (c_mid - u_mid**2) - m3 * v_mid
            R4 = dg_deta - p_mid
            R5 = dep_deta + dduv_deta + m1 * f_mid * p_mid  - m3 * p_mid

            #Columns of A matrix
            i0 = 5 * j    # Corresponds to j
            i1 = 5 * (j+1)          # Corresponds to j+1

            # R1: (f[j+1]-f[j])/hj - u_jmid = 0
            row = 3 + 5*j
            A_mat[row, i0+0] = -1.0/hj
            A_mat[row, i0+1] = -0.5
            A_mat[row, i1+0] = +1.0/hj
            A_mat[row, i1+1] = -0.5
            b_vec[row] = -R1

            # R2: (u[j+1]-u[j])/hj - v_jmid = 0
            row += 1
            A_mat[row, i0+1] = -1.0/hj
            A_mat[row, i0+2] = -0.5
            A_mat[row, i1+1] = +1.0/hj
            A_mat[row, i1+2] = -0.5
            b_vec[row] = -R2

            # R3 partial derivatives
            dR3_df_j    = 0.5 * m1 * v_mid
            dR3_df_jp1  = 0.5 * m1 * v_mid
            dR3_du_j    = - m2 * u_mid
            dR3_du_jp1  = - m2 * u_mid
            dR3_dv_j    = - b[j]/hj + 0.5 * m1 * f_mid - 0.5 * m3
            dR3_dv_jp1  = b[j+1]/hj + 0.5 * m1 * f_mid - 0.5 * m3

            row += 1
            A_mat[row, i0+0] = dR3_df_j
            A_mat[row, i0+1] = dR3_du_j
            A_mat[row, i0+2] = dR3_dv_j
            A_mat[row, i1+0] = dR3_df_jp1
            A_mat[row, i1+1] = dR3_du_jp1
            A_mat[row, i1+2] = dR3_dv_jp1
            b_vec[row] = -R3

            # R4 (g[j+1]-g[j])/hj - p_jmid = 0
            row += 1
            A_mat[row, i0+3] = -1.0/hj
            A_mat[row, i0+4] = -0.5
            A_mat[row, i1+3] = +1.0/hj
            A_mat[row, i1+4] = -0.5
            b_vec[row] = -R4

            #R5 partial derivatives
            dR5_df_j    = 0.5 * m1 * p_mid 
            dR5_df_jp1  = 0.5 * m1 * p_mid 
            dR5_du_j    = - (d[j] * v[j]) / hj 
            dR5_du_jp1  = (d[j+1] * v[j+1]) / hj 
            dR5_dv_j    = - (d[j] * u[j]) / hj
            dR5_dv_jp1  = (d[j+1] * u[j+1]) / hj
            dR5_dg_j    = 0
            dR5_dg_jp1  = 0
            dR5_dp_j    = - e[j] / hj + 0.5 * m1 * f_mid - 0.5 * m3 
            dR5_dp_jp1  = e[j+1] / hj + 0.5 * m1 * f_mid - 0.5 * m3 

            row += 1
            A_mat[row, i0+0] = dR5_df_j
            A_mat[row, i0+1] = dR5_du_j
            A_mat[row, i0+2] = dR5_dv_j
            A_mat[row, i0+3] = dR5_dg_j
            A_mat[row, i0+4] = dR5_dp_j
            A_mat[row, i1+0] = dR5_df_jp1
            A_mat[row, i1+1] = dR5_du_jp1
            A_mat[row, i1+2] = dR5_dv_jp1
            A_mat[row, i1+3] = dR5_dg_jp1
            A_mat[row, i1+4] = dR5_dp_jp1
            b_vec[row] = -R5


        #Boundary conditions at freestream u(eta) = 1 and g(eta) = 1 
        A_mat[-2, 5*(J-1)+1] = 1.0
        b_vec[-2] = 1.0 - u[J-1]

        A_mat[-1, 5*(J-1)+3] = 1.0
        b_vec[-1] = 1.0 - g[J-1]

        delta = np.linalg.solve(A_mat, b_vec)

        for j in range(J):
            f[j] += delta[5 * j + 0]
            u[j] += delta[5 * j + 1]
            v[j] += delta[5 * j + 2]
            g[j] += delta[5 * j + 3]
            p[j] += delta[5 * j + 4]

        res_norm = max(abs(R1), abs(R2), abs(R3), abs(R4), abs(R5))       
        if np.max(np.abs(delta)) < 1e-6 and res_norm < 1e-6:
            print(f"Solution at x=0 converged after {it} iterations")
            break
    
    plt.figure()
    plt.plot(u, eta)
    plt.show() 
    
    bl = BL(f=f, u=u, v=v, g=g, p=p)
    return bl



def march_downstream(grid, bl, pressure, max_iter):
    xi = grid.xi
    k = grid.k
    eta = grid.eta
    h = grid.h
    N = len(xi)
    J = len(eta)

    #Use the initial f, u, v as xi-1 for the first iteration
    f_nm1 = np.copy(bl.f)
    u_nm1 = np.copy(bl.u)
    v_nm1 = np.copy(bl.v)
    g_nm1 = np.copy(bl.g)
    p_nm1 = np.copy(bl.p)

    #Pressure Gradient arrays
    m1_array = pressure.m1
    m2_array = pressure.m2
    m3_array = pressure.m3

    for n in range(1, N):

        kn = k[n-1]
        xi_mid = xi[n] - kn/2

        #Use values of xi point before as initial guess
        f = np.copy(f_nm1)
        u = np.copy(u_nm1)
        v = np.copy(v_nm1)
        g = np.copy(g_nm1)
        p = np.copy(p_nm1)

        m1 = m1_array[n]
        m2 = m2_array[n]
        m3 = m3_array[n]

        for it in range(1, max_iter):

            #Create matrix A and vector b for system of equations
            size  = 5 * J
            A_mat = np.zeros((size, size))
            b_vec = np.zeros(size)

            #First two rows for boundary conditions f(0) = 0 and u(0) = 0
            A_mat[0, 0] = 1.0
            b_vec[0]    = 0.0 - f[0]

            A_mat[1, 1] = 1.0
            b_vec[1]    = 0.0 - u[0]

            #g(0) = gw
            A_mat[2, 3] = 1.0
            b_vec[2] = g_w - g[0]

            b, c, d, e = compute_local_coefficients(eta, f, u, v, g, p, xi=xi_mid)

            for j in range(0, J-1):
                hj = h[j]

                #Eta Derivatives
                df_deta = (f[j+1] - f[j])/hj
                du_deta = (u[j+1] - u[j])/hj
                dg_deta = (g[j+1] - g[j])/hj
                dbv_deta = (b[j+1]*v[j+1] - b[j]*v[j])/hj
                dep_deta = (e[j+1]*p[j+1] - e[j]*p[j])/hj
                dduv_deta = (d[j+1]*u[j+1]*v[j+1] - d[j]*u[j]*v[j])/hj

                #Mid Points (j+1)
                f_mid = (f[j] + f[j+1])/2
                f_nm1_mid = (f_nm1[j] + f_nm1[j+1])/2
                u_mid = (u[j] + u[j+1])/2
                u_nm1_mid = (u_nm1[j] + u_nm1[j+1])/2
                v_mid = (v[j] + v[j+1])/2
                b_mid = (b[j] + b[j+1])/2
                c_mid = (c[j] + c[j+1])/2
                p_mid = (p[j] + p[j+1])/2
                g_mid = (g[j] + g[j+1])/2
                g_nm1_mid = (g_nm1[j] + g_nm1[j+1])/2

                #Xi Derivatives
                df_mid_dxi = (f_mid - f_nm1_mid)/kn
                du_mid_dxi = (u_mid - u_nm1_mid)/kn
                dg_mid_dxi = (g_mid - g_nm1_mid)/kn

                #Residuals
                R1 = df_deta - u_mid
                R2 = du_deta - v_mid
                R3 = dbv_deta + m1 * f_mid * v_mid + m2 * (c_mid - u_mid**2) - m3 * v_mid - xi_mid * (u_mid * du_mid_dxi - v_mid * df_mid_dxi)
                R4 = dg_deta - p_mid
                R5 = dep_deta + dduv_deta + m1 * f_mid * p_mid  - m3 * p_mid - xi_mid * (u_mid * dg_mid_dxi - p_mid * df_mid_dxi)

                #Columns of A matrix
                i0 = 5 * j    # Corresponds to j
                i1 = 5 * (j+1)          # Corresponds to j+1

                # R1: (f[j+1]-f[j])/hj - u_jmid = 0
                row = 3 + 5*j
                A_mat[row, i0+0] = -1.0/hj
                A_mat[row, i0+1] = -0.5
                A_mat[row, i1+0] = +1.0/hj
                A_mat[row, i1+1] = -0.5
                b_vec[row] = -R1

                # R2: (u[j+1]-u[j])/hj - v_jmid = 0
                row += 1
                A_mat[row, i0+1] = -1.0/hj
                A_mat[row, i0+2] = -0.5
                A_mat[row, i1+1] = +1.0/hj
                A_mat[row, i1+2] = -0.5
                b_vec[row] = -R2

                # R3 partial derivatives
                dR3_df_j    = 0.5 * m1 * v_mid + 0.5 * xi_mid * v_mid / kn
                dR3_df_jp1  = 0.5 * m1 * v_mid + 0.5 * xi_mid * v_mid / kn
                dR3_du_j    = - m2 * u_mid - xi_mid * (0.5 * du_mid_dxi + 0.5 * u_mid / kn )
                dR3_du_jp1  = - m2 * u_mid - xi_mid * (0.5 * du_mid_dxi + 0.5 * u_mid / kn )
                dR3_dv_j    = - b[j]/hj + 0.5 * m1 * f_mid - 0.5 * m3 + 0.5 * xi_mid * df_mid_dxi
                dR3_dv_jp1  = b[j+1]/hj + 0.5 * m1 * f_mid - 0.5 * m3 + 0.5 * xi_mid * df_mid_dxi

                row += 1
                A_mat[row, i0+0] = dR3_df_j
                A_mat[row, i0+1] = dR3_du_j
                A_mat[row, i0+2] = dR3_dv_j
                A_mat[row, i1+0] = dR3_df_jp1
                A_mat[row, i1+1] = dR3_du_jp1
                A_mat[row, i1+2] = dR3_dv_jp1
                b_vec[row] = -R3

                # R4 (g[j+1]-g[j])/hj - p_jmid = 0
                row += 1
                A_mat[row, i0+3] = -1.0/hj
                A_mat[row, i0+4] = -0.5
                A_mat[row, i1+3] = +1.0/hj
                A_mat[row, i1+4] = -0.5
                b_vec[row] = -R4

                #R5 partial derivatives
                dR5_df_j    = 0.5 * m1 * p_mid + 0.5 * xi_mid * p_mid /kn
                dR5_df_jp1  = 0.5 * m1 * p_mid + 0.5 * xi_mid * p_mid /kn
                dR5_du_j    = - (d[j] * v[j]) / hj - 0.5 * xi_mid * dg_mid_dxi
                dR5_du_jp1  = (d[j+1] * v[j+1]) / hj  - 0.5 * xi_mid * dg_mid_dxi
                dR5_dv_j    = - (d[j] * u[j]) / hj
                dR5_dv_jp1  = (d[j+1] * u[j+1]) / hj
                dR5_dg_j    = - 0.5 * xi_mid * u_mid / kn
                dR5_dg_jp1  = - 0.5 * xi_mid * u_mid / kn
                dR5_dp_j    = - e[j] / hj + 0.5 * m1 * f_mid - 0.5 * m3 + 0.5 * xi_mid * df_mid_dxi
                dR5_dp_jp1  = e[j+1] / hj + 0.5 * m1 * f_mid - 0.5 * m3 + 0.5 * xi_mid * df_mid_dxi

                row += 1
                A_mat[row, i0+0] = dR5_df_j
                A_mat[row, i0+1] = dR5_du_j
                A_mat[row, i0+2] = dR5_dv_j
                A_mat[row, i0+3] = dR5_dg_j
                A_mat[row, i0+4] = dR5_dp_j
                A_mat[row, i1+0] = dR5_df_jp1
                A_mat[row, i1+1] = dR5_du_jp1
                A_mat[row, i1+2] = dR5_dv_jp1
                A_mat[row, i1+3] = dR5_dg_jp1
                A_mat[row, i1+4] = dR5_dp_jp1
                b_vec[row] = -R5


            #Boundary conditions at freestream u(eta) = 1 and g(eta) = 1 
            A_mat[-2, 5*(J-1)+1] = 1.0
            b_vec[-2] = 1.0 - u[J-1]

            A_mat[-1, 5*(J-1)+3] = 1.0
            b_vec[-1] = 1.0 - g[J-1]

            delta = np.linalg.solve(A_mat, b_vec)

            for j in range(J):
                f[j] += delta[5 * j + 0]
                u[j] += delta[5 * j + 1]
                v[j] += delta[5 * j + 2]
                g[j] += delta[5 * j + 3]
                p[j] += delta[5 * j + 4]
            
            res_norm = max(abs(R1), abs(R2), abs(R3), abs(R4), abs(R5))
            if np.max(np.abs(delta)) < 1e-6 and res_norm < 1e-6:
                print(f"Solution at x={xi[n]} converged after {it} iterations")
                break

        f_nm1 = np.copy(f)
        u_nm1 = np.copy(u)
        v_nm1 = np.copy(v)




def main():
    max_iter = get_input()
    
    grid = build_grid(xi_max=10, NX = 11, eta_max=8, NY=501)

    pressure_coeff = calc_pressure_coefficients(grid)

    bl = initial_bl(grid, pressure_coeff, max_iter)

    bl_thickness = march_downstream(grid, bl, pressure_coeff, max_iter)


main()