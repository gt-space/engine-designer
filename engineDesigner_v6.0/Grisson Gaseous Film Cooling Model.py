import numpy as np
import scipy as sp

#gaseous grisson model

'''
Mbl = mass flux of the boundary layer
Taw = adiabatic wall temperature
Tr = recovery temperature
Cpg = specific heat at constant pressure of the gas
Cvg = specific heat at constant volume of the gas
Cpc = specific heat at constant pressure of the coolant
mWg = molecular weight of the gas
mWc = molecular weight of the coolant
Pr = Prandtl number
Tc = temperature of the coolant
Tu = turbulence coefficient
Ts = static temperature of the gas
Gvs = gas viscosity
x = distance along the surface **FROM THE DRYOUT POINT**
r = local radius
density = core gas density
v = core gas velocity
St = Stanton number
'''

gamma = Cpg / Cvg

n = len(x)
diameter = 2*r
gFlux = density * v

#solving for local mach number of the free stream gas
gasConstant = Cpg * (gamma - 1) / gamma
mach = v / np.sqrt(gamma * gasConstant * Ts)

#recovery temperature is the adiabatic wall temp 
recoveryFactor = float(Pr ** (1/3))
Tr = Ts * (1 + recoveryFactor * (gamma - 1) / 2 * mach ** 2)

#this a grisson empirical factor where Tu = 0.1 
ampFactor = 1 + 10.2 * Tu
#the amp factor tells us how much the heat transfer is amplified cause of turbulence

#foreign gas factor
if mWg > mWc:
    fgf = (mWc / mWg) ** 0.60
else:
    #a heavier coolant has a weaker effect
    fgf = (mWg / mWc) ** 0.35
#the foreign gas factor "fixes" the difference in molecular weight between the gas and the coolant

Mbl = np.zeros(n)
Taw = np.zeros(n)

conservedMass = np.zeros(n)

#the initial mass flow at the axial injection point
baseMdot = coolantMdot / (np.pi * diamater[0])

Mbl[0] = baseMdot
Taw[0] = Tc
conservedMass[0] = baseMdot * diameter[0]

def ODE(y, massConversedVal, TawVal):
    #need to interpolate to get the values at a certain y (axial distance)
    diameterInt = float(np.interp(y, x, diameter))
    gFluxInt = float(np.interp(y, x, gFlux))
    TrInt = float(np.interp(y, x, Tr))
    gvsInt = float(np.interp(y, x, Gvs))

    #this is the Mbl at the current y value using the conserved mass value
    mbl = massConversedVal / diameterInt

    #the 0.1963 is a grisson empirical constant for turbulence skin friction
    #this equation is the entrainment relationship grisson uses to calculate the change in mass flow of the boundary layer with respect to x
    dMbldx = 0.1963 * ampFactor * gFluxInt * (gvsInt / mbl) ** 0.25
    dMassConverseddx = dMbldx * diameterInt

    #this next equation, the tempeature ODE, is split into two parts, the numerator and demoninator
    # the numerator is the change in temperature of the adiabatic wall with respect to x
    # the denominator is the effective thermal mass of the boundary layer

    denominator = mbl + (Cpc / (fgf * Cpg) - 1) * baseMdot
    numerator = dMbldx * (TrInt - TawVal)

    dTawdx = numerator / denominator


    return dTawdx, dMassConverseddx

#this is a slight deviation from the original grisson model
# instead of the original model which prefers the finite difference method to solve the ODEs
# I'm using the RK4 approximation method
for i in range(n - 1):
    xi = x[i]
    dx = x[i+1] - x[i]

    cMass = conservedMass[i]
    taw = Taw[i]

    #now must get the k values for the RK4 method
    k1Taw, k1Mass = ODE(xi, cMass, taw)

    midpointX = xi + dx / 2
    midpointcMass = cMass + k1Mass * dx / 2
    midpointTaw = taw + k1Taw * dx / 2
    k2Taw, k2Mass = ODE(midpointX, midpointcMass, midpointTaw)

    midpointcMass2 = cMass + k2Mass * dx / 2
    midpointTaw2 = taw + k2Taw * dx / 2
    k3Taw, k3Mass = ODE(midpointX, midpointcMass2, midpointTaw2)

    endpointX = xi + dx
    endpointcMass = cMass + k3Mass * dx
    endpointTaw = taw + k3Taw * dx
    k4Taw, k4Mass = ODE(endpointX, endpointcMass, endpointTaw)

    #now we can calculate the next values for mass and taw
    conservedMass[i+1] = cMass + (dx / 6) * (k1Mass + 2 * k2Mass + 2 * k3Mass + k4Mass)
    Taw[i+1] = taw + (dx / 6) * (k1Taw + 2 * k2Taw + 2 * k3Taw + k4Taw)

    Mbl[i+1] = conservedMass[i+1] / diameter[i+1]

effectiveness = (Tr - Taw) / (Tr - Tc)
effectiveness[effectiveness < 0] = 0
effectiveness[effectiveness > 1] = 1

#this is a turbulent velocity profile used in grisson's model for the boundary layer mass flux
# here, we are rearranging the equation to figure out the ending boundary layer thickness
boundaryLayerThickness = (8/7) * (Mbl / gFlux)
endingBoundaryLayerThickness = changeInMbl[-1]




    

