import scipy as sp
import numpy as np

'''
Free stream mass flux = G
specific heat gas = Cpg
specific heat coolant = Cpc
chamber diameter = D
gas viscosity = Gvs
prandtl number = Pr
stag temperature = To
recovery temp = Tr
adibatic wall temp = Tad
coolant saturation temp = Tsat
initial coolant temp = Tcf
latent heat of vaporization = hfg
coolant mass flow rate = m_dot_coolant
liquid film length = fcl
rp1 specific heat at constant pressure = Cp
rp1 specific heat at constant volume = Cv
mach number = M

this stuff isn't necessarily needed but it is useful for calculating the heat fluxes
boltzmann constant = *scipy.constants.sigma*
effective emissivity = gasEp
absorptivity of film = lAbs
absorptivity of wall = wAbs
heat flux limit = burnoutLimit
'''
gamma = Cp / Cv

#static temperature
def staticTemperature(G, D, Gvs, Cpg, Pr, M, gamma):
    staticTemp = To / (1 + ((gamma - 1) / 2) * M**2)
    return staticTemp

#recovery paramter + temperature
def recoveryParameter(Pr):
    r = Pr**(1/3)
    return r
def recoveryTemperature(staticTemp, r, gamma, M):
    Tr = staticTemp * (1 + r * ((gamma - 1) / 2) * M**2)
    return Tr

#bartz w/ stanton number correlation
def bartz(G, D, Gvs, Cpg, Pr):
    h0 = 0.026 * (G**0.8) * (D**-0.2) * (Gvs**-0.2) * Cpg * (Pr**-0.6)
    return h0

#blowing parameter
def blowingParamter(Cpg, Tr, Tsat, hfg, Cpc, Tcf):
    B = (Cpg * (Tr - Tsat)) / (hfg + Cpc * (Tsat - Tcf))
    return B

def adjustedh0(h0, B):
    hg = h0 / (np.log(1 + B) / B)
    return hg

constants = dict(r=recoveryParameter(Pr),
                 staticTemp=staticTemperature(G, D, Gvs, Cpg, Pr, M, gamma),
                 recoveryTemp=recoveryTemperature(staticTemp, r, gamma, M),
                 h0=bartz(G, D, Gvs, Cpg, Pr),
                 B=blowingParamter(Cpg, Tr, Tsat, hfg, Cpc, Tcf),
                 hg=adjustedh0(h0, B))

hfgActual = hfg + Cpc * (Tsat - Tcf)


#grisson liquid ffc

#convection heat fluxes
filmHeatFluxConv = hg * (Tr - Tsat)
wallHeatFluxConv = h0 * (Tr - Tad)

#radiation heat fluxes
filmHeatFluxRad = gasEp * sp.constants.sigma * (To**4 - Tsat**4) * lAbs
wallHeatFluxRad = gasEp * sp.constants.sigma * (To**4 - Tad**4) * wAbs

#actual heat fluxes
filmHeatFlux = filmHeatFluxConv + filmHeatFluxRad
wallHeatFlux = wallHeatFluxConv + wallHeatFluxRad

#step count
n = 1001

#step size + x array
dx = chamberLength / (n - 1)
x_arr = np.linspace(0, chamberLength, n)


#initialized arrays
m_dot_coolant_arr = np.zeros(n)
Tfilm = np.zeros(n)
qFilm = np.zeros(n)
qWall = np.full(n, wallHeatFlux)

dryoutPoint = 0
burnoutPoint = 0
burnedOut = False


chamberRadius = D / 2
tempor_m_dot_coolant = m_dot_coolant

for i in range(n):
    x = x_arr[i]

    if tempor_m_dot_coolant > 0 and not burnedOut:
        if filmHeatFlux > burnoutLimit:
            burnedOut = True
            burnoutPoint = x

    valid = tempor_m_dot_coolant > 0 and not burnedOut

    if valid: #if the film is there
        q = filmHeatFlux
        temp = Tsat
    else: #if the film has evaporated
        q = wallHeatFlux
        temp = Tr 
        if dryoutPoint == 0 and not burnedOut: #if the film just evaporated
            dryoutPoint = x

    m_dot_coolant_arr[i] = max(tempor_m_dot_coolant, 0)
    Tfilm[i] = temp
    qFilm[i] = q

    if valid:
        dm_mdot = (-( 2 * np.pi * chamberRadius * q ) / hfgActual) * dx
        tempor_m_dot_coolant += dm_mdot

if dryoutPoint == 0: #film somehow makes it the the throat without evaporating
    dryoutPoint = chamberLength

effiency_arr = (Tr - Tfilm) / (Tr - Tcf)
meanEfficiency = np.mean(effiency_arr)
#liquidEfficiency = (Tad - Tadf) / (Tad - Tc, 0)
