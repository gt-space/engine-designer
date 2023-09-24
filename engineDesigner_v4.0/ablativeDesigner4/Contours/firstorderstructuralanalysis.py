# Inputs
FOS = 3
DT_OW = 53
p_OW = 400 * 6894.76 #Pa
r_OW = 2.02 * 0.0254 #m
# Mat props: OW
k_OW = 10 # conservative thing https://www.researchgate.net/profile/Ronald-Joven/publication/288102626_Thermal_properties_of_carbon_fiber-epoxy_composites_with_different_fabric_weaves/links/56be27e408aee5caccf2f5d3/Thermal-properties-of-carbon-fiber-epoxy-composites-with-different-fabric-weaves.pdf
rho_OW = 1410  # kg/m3
Cp_OW = 1200 # J/kg*K
v_OW = .4  # using epoxy value
E_OW = 7e10  # Pa Matweb
CTE_OW = 2.13e-5 #m/m/C = m/m/K
v = .4
E_OW = 7e10
UTS_OW = 918e6 #Pa Matweb
# Fxns for OverWrap
Sig_th = 2 * CTE_OW * E_OW * DT_OW / (1 - v)
t_OW = p_OW * r_OW / (UTS_OW / FOS - Sig_th)

# Mat props: Insulator
k_I = .225  # W/m*K #number from literature
rho_I = 2300  #
Cp_I = 1000  #
CTE_I = 17e-6  #
v_I = .23  #
E_I = 36e9  #
UTS_I = 15.2e6  #

# Insulator FXNs
t_I_perDT_Ins = (k_I/k_OW)*t_OW/DT_OW
t_I_max = .0625*.0254
max_DT_Ins= t_I_max/t_I_perDT_Ins

#q_dot_solver
# let
T_hg = 5000

def bartz(engine, T_wg, i):
    # DUE TO BARTZ CORRELATION BEING DEVELOPED IN ENGLISH UNITS, VALUES
    # CONVERTED TO ENGLISH UNITS IN EQUATIONS BELOW THEN RECONVERTED
    # THIS IS THE BARTZ

    recovery_factor = 1

    g = 9.8104 # Gravity, m/s^2
    C_star = engine.C_star
    R_tCurve = engine.R_tCurve
    R_t = engine.R_t
    P_c = engine.P_inj_psi * 6894.76 # Injector face pressure in psi
    T_c0 = engine.engineProps[0,9] # Convert K to R
    gam0 = engine.engineProps[0,15] # Chamber stagnation gamma
    # Average of frozen and not frozen
    cp_ns = (engine.engineProps[99,20] + engine.engineProps[99,14]) * 0.5 * 1000#Convert from Joules to BTU/lb*F
    praneff_ns = 0.5 * (engine.engineProps[99,22] + engine.engineProps[99,19]) # No conversion needed
    '''
    praneff_ns_f = engine.engineProps[99,22] # No conversion needed
    praneff_ns_eq = engine.engineProps[99,19]

    visc_ns = (engine.engineProps[0,17]/1000) * 0.1 # Dyn viscosity (mili-Poise > Poise > kg/ms)
    cond = 0.1 * (engine.engineProps[i,21] + engine.engineProps[i,18]) * 0.5 # Eq Conductivity (W/mK)
    cp_ns_calc = (engine.engineProps[99,20] + engine.engineProps[99,14]) * 1000 #(kJ > J)
    # praneff_ns = visc_ns*cp_ns_calc/cond
    # praneff_ns = 0.508
    # praneff_ns_est = (4*gam0)/(9*gam0-5) #Originally used effective prandtl number from CEA, switched to gamma correlation
    '''
    visc_ns = (engine.engineProps[0,17]/1000) * 0.1# Convert to Poise

    # Run Bartz Correlation
    gam = engine.engineProps[i, 15] #Index gamma; No conversion needed
    mach = engine.engineProps[i, 4]
    contourR = engine.engineProps[i, 0] #Convert to Inches
    T_aw = T_c0 * (1 + (praneff_ns**(1/3)) * ((gam - 1)/2) * (mach**2))/(1 + ((gam - 1)/2) * (mach**2))
    T_aw = T_aw * recovery_factor
    #Sigma = sigA * sigB - Split into Terms for readability
    sigA = (0.5 * (T_wg/T_c0) * (1 + ((gam - 1)/2) * mach**2) + 0.5)**(-0.68)
    sigB = (1 + ((gam - 1)/2) * (mach**2))**(-0.12)
    sigma = sigA*sigB
    # a_T = (aA*aB*aC*aD)*aE*sigma - Split into Terms
    aA = 0.026/((2 * R_t)**0.2)
    aB = ((visc_ns**0.2) * cp_ns)/(praneff_ns**0.6)
    aC = (P_c/C_star)**0.8
    aD = ((2 * R_t)/R_tCurve)**0.1
    aE = ((R_t**2)/(contourR**2))**0.9
    h_g = (aA * aB * aC * aD) * aE * sigma #Convection coefficient
    q_conv_max = h_g * (T_aw - T_wg) #Heat flux (heat per area)

    # Adiabtic wall temp T_aw is temp if the wall was adiabatic (no heat transfer
    # this means that this would be equal to the flow temperature at that region
    # thus, T_aw = T_hg

    # print(h_g)
    return (h_g, q_conv, T_aw)


(h_g, q_conv, T_aw) = bartz(self.engine, T_w, wall_idx)

def q(self,T_w,T_hg,wall_idx):
    return bartz(self.engine,T_w, wall_idx)* (T_hg - T_w)


def dq_dT_w(T_w, T_hg,wall_idx):
    # Calculate the derivative of q with respect to T_w
    epsilon = 1e-6  # Small value for numerical differentiation
    dq = self.q(T_w,T_hg,wall_idx)
    dq_dT_w = (self.q(T_w + epsilon, T_hg,wall_idx) - dq) / epsilon
    return dq_dT_w


def newton_raphson(q_target, T_hg, initial_guess, wall_idx, max_iterations=100, tolerance=1e-6):
    T_w = initial_guess
    for iteration in range(max_iterations):
        q_current = self.q(T_w, T_hg,wall_idx)
        dq_dT_w_current = self.dq_dT_w(T_w, T_hg,wall_idx)

        T_w_next = T_w - (q_current - q_target) / dq_dT_w_current

        if abs(T_w_next - T_w) < tolerance:
            return T_w_next

        T_w = T_w_next

    raise RuntimeError("Newton-Raphson did not converge.")


# Provide your known values for q_target and T_hg
q_target = ...
T_hg = ...

# Provide an initial guess for T_w
initial_guess = ...

# Call the newton_raphson function to find the solution for T_w
solution = newton_raphson(q_target, T_hg, initial_guess)

print("Solution for T_w:", solution)

#Goal: Find T_w that gives q_dot that aligns with the prescribed DT based on structural overwrap
#do newton raphson to do this

#then solve for qdot

# then use qdot to find TII and TIO

# Use this to find t_I


# Outputs
print("Overwrap Thickness [in]: " + str(t_OW/.0254))
print("Max Insulator DT given max thickness of "+str(t_I_max/.0254) +" in: " + str(max_DT_Ins)+str(" K"))