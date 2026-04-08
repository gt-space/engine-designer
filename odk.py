#----------------------------------------------------------------
#           ONE-DIMENSIONAL KINETICS (ODK) MODULE
# --------------------------------------------------------------- 

# I haven't decided how the class is going to work yet 

import numpy as np

#things we're going to have to input somehow

V = []
dCidx = []
gamma = []
dgammadT = []
dgammadCi = []

Rd = []
x = []
Y = []
dYdx = []
d2Ydx2 = []

wi = []
rStar = []
rho = []

R = []
Ri = []

T = []
Cpi = []

hi = []


phi1 = []
phi2 = []
phi3 = []

P = []
dPdx = []
d2Pdx2 = []
  
"""The NASA report has functionality for conical nozzles as well,
   but we don't really have a reason to add it."""

"""In this equation, the concentration of a species is Ci."""
dCidx = wi * rStar / (rho * V)
  
S1 = 1/R * sum(dCidx*Ri)
dS1dV = 1/R * sum(phi1*Ri)
dS1drho = 1/R * sum(phi2*Ri)
dS1dT = 1/R * sum(phi3 * Ri)
dS1dCi = [...]
  
S2 = 1/R/T * sum(dCidx*hi)
dS2dV = 1/R/T * sum(phi1*hi)
dS2drho = 1/R/T * sum(phi2 * hi)
dS2dT = 1/R/T * sum(phi3 * hi + dCidx*Cpi) - S2/T
dS2dCi = [...]
   

gRat = (gamma-1)/gamma
B = gRat * S2
dBdV = gRat * dS2dV
dBdrho = gRat * dS2drho
dBdT = gRat * dS2dT + S2 * dgammadT/gamma^2
dBdCi = gRat * dS2dCi + S2 * dgammadCi/gamma^2

A = S1- B
dAdV = dS1dV -  dBdV
dAdrho = dS1drho - dBdrho
dAdT = dS1dT - dBdT
dAdCi = dS1dCi - dBdCi

Msq = V^2/gamma/R/T
dMsqdV = 2*Msq/V
dMsqdT = -Msq/T - Msq/gamma * dgammadT
dMsqdCi = -Msq * (dgammadCi/gamma + Ri/R)


def subBetaDict(rho, V, dPdx, d2Pdx2, P, A, gamma, T, B):
   
   gRat = (gamma-1)/gamma

   dVdx = -1/rho/V * dPdx
   d2Vdx2 = -1/rho/V * d2Pdx2
   drhodx = rho * (dPdx * 1/P/gamma - A)
   d2rhodx2 = rho/gamma/P * (d2Pdx2 - (dPdx)**2 * 1/P)
   dTdx = T * (gRat/P*dPdx - B)
   d2Tdx2 = gRat * T/P * (d2Pdx2-(dPdx)**2 / P)

   subBeta = [{"V":-1/V * dVdx, "rho":-1/rho * dVdx}]
   subBeta.append({"V":-rho * dAdV, "rho": -1/rho * drhodx - rho*dAdrho, "T":-rho*dAdT-rho/P/gamma**2*dgammadT*dPdx, "Ci":-rho/gamma**2/P*dgammadCi*dPdx-rho*dAdCi})
   subBeta.append({"V":-T*dBdV, "rho": -T*dBdrho, "T": 1/T*dTdx + T/gamma**2/P*dPdx*dgammadT-T*dBdT, "Ci": T * (1/gamma**2/P*dPdx*dgammadCi-dBdCi)})
   
   return subBeta

def aRats(Rd, x):
   
   at = (1+Rd-np.sqrt(Rd^2-x^2))^2
   datdx = 2*x*(1+Rd-np.sqrt(Rd^2-x^2))/np.sqrt(Rd^2-x^2)
   d2atdx2 = (2/np.sqrt(Rd^2-x^2) + 2*x^2/(Rd-x^2)^1.5) * (1+Rd-np.sqrt(Rd-x^2)) + 2*x^2/Rd^2-x^2
   
   ac = Y^2
   dacdx = 2 * Y * dYdx
   d2acdx2 = 2 * (Y*d2Ydx2 + dYdx^2)
   
   return 


def supBetaDict(a, dadx, d2adx2, V, M):
   dVdx = V/(Msq-1) * (1/a * dadx - A)
   d2Vdx2 = V/(Msq-1) * 1/a * (d2adx2 - 1/a * dadx**2)
   
   drhodx = -rho * (Msq/(Msq-1) * (dadx/a-A)+A)
   d2rhodx2 = -rho * Msq/(Msq-1) * 1/a * (d2adx2-1/a * (dadx)**2)

   dTdx = -T * ((gamma-1) * Msq/(Msq-1)*(dadx/a-A)+B)
   d2Tdx2 = -T*Msq/(Msq-1) * (gamma-1)/a * (d2adx2-1/a*dadx**2)

   supBeta = [{"V": 1/V*dVdx-1/(Msq-1)*dVdx*dMsqdV-V/(Msq-1)*dAdV, "rho": -V/(Msq-1)*dAdrho, "T": -1/(Msq-1) * dVdx * dMsqdT - V/(Msq-1) * dAdT, "Ci": -1/(Msq-1)*dVdx*dMsqdCi-V/(Msq-1)*dAdCi}]
   supBeta.append({"V": rho*((1/(Msq-1)**2)*(1/a*dadx-A)*dMsqdV+dAdV/(Msq-1)), "rho": 1/rho*drhodx + rho/(Msq-1)*dAdrho, "T": rho * ((1/(Msq-1)**2)*(dadx/a-A)*dMsqdT+dAdT/(Msq-1)), "Ci": rho * ((dadx/a/-A)*dMsqdCi/(Msq-1)**2+dAdCi/(Msq-1))})
   supBeta.append({"V": T * ((gamma-1)/(Msq-1)**2 * (dadx/a-A) * dMsqdV + (gamma-1)*Msq/(Msq-1)*dAdV-dBdV), "rho": T * ((gamma-1)*Msq/(Msq-1)*dAdrho-dBdrho), "T": dTdx/T + T*((gamma-1)/(Msq-1)**2*(dadx/a-A)*dMsqdT + (gamma-1)*Msq/(Msq-1)*dAdT-dBdT), "Ci": T * ((gamma-1)/(Msq-1)**2*(dadx/a-A)*dMsqdCi + (gamma-1)*Msq/(Msq-1)*dAdCi-dBdCi-Msq/(Msq-1)*(dadx/a-A)*dgammadCi)})

   return supBeta



# need to work on iteration and then we're good to go