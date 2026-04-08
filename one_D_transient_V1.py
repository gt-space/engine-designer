import cantera as ct
import math
import matplotlib.pyplot as plt
import numpy as np
import time
from pathlib import Path
script_dir = Path(__file__).parent

#conversions
psiaToPa = 6894.76
inch_to_meter = 0.0254

#initial conditions
mfuel = 0.08
mlox = 0.6151
mtotal = mfuel + mlox
chamber_length = 9 * inch_to_meter
chamber_radius = 2.875 * inch_to_meter
area = (chamber_radius**2)*math.pi
steps = 100 #increase number of steps for both soon
timeSteps = 10
dt = 1 #make smaller later when number of timeSteps is increased
dz = chamber_length/steps
dV = dz*area

#reactor 1 setup
initialMixture = ct.Solution("nDodecane_Reitz.yaml")
initialMixture.TP = 1500, 300*psiaToPa #sets 300 psia as expected initial chamber pressure; temperature is a placeholder --> by the end, if this works right, it should reach equilibrium (given enough time steps?)
initialMixture.set_mixture_fraction(0.338,"C12H26:1","O2:1")

reactor1 = ct.Reactor(initialMixture, clone=True)
reactor1.volume = dV

#reservoirs representing the reactor directly upstream/downstream --> upstream refers to the lox/fuel inlets, and downstream reservoir is technically irrelevant (same as mixture))

#lox inlet
#lox = ct.Solution('nDodecane_Reitz.yaml')
#lox.TPX = 291.758, 363.594*psiaToPa, 'o2:1' 
#res_lox = ct.Reservoir(lox, name="LOX Reservoir", clone=True) 

#fuel inlet
#fuel = ct.Solution('nDodecane_Reitz.yaml')
#fuel.TPX = 400, 751.906*psiaToPa, 'c12h26:1.0'
#res_fuel = ct.Reservoir(fuel, name="Fuel Reservoir", clone=True) 
upstream = ct.Reservoir(initialMixture, name="upstream", clone=True)
downstream = ct.Reservoir(initialMixture, name="downstream", clone=True) #what reactor 1 will exhaust into

#flow into reactor 1 --> fuel and lox inlets
#FMFC = ct.MassFlowController(res_fuel, reactor1, mdot = mfuel, name = "Fuel Inlet") #from previous (0-dimensional reactors) file
#OMFC = ct.MassFlowController(res_lox, reactor1, mdot = mlox, name = "LOX inlet") #idk if these are necessary yet

m = ct.MassFlowController(upstream, reactor1, mdot=mtotal)

#flow out of reactor 1 AKA flow into reactor 2 --> this means it also establishes dP
v = ct.PressureController(reactor1, downstream, primary=m, K=1e-12) #value of K affects transient flow but idk how to establish it
#note on "primary" --> if the pressure controller can take two values for upstream (FMFC and OMFC), i think it might be better to switch to that --> for now, i'm sticking to how the cantera reference set this up

#set up network
network = ct.ReactorNet([reactor1])
network.max_time_step = 1e4 #same value from PFR page on cantera --> seems like a very high limit tho? --> i don't use this anymore nvm


#setup for loop 
z = (np.arange(steps) + 1) * dz
t_reactor1 = np.zeros_like(z)  # residence time in each reactor
u = np.zeros_like(z)
t = np.zeros_like(z)
states = ct.SolutionArray(initialMixture)

#plotting 
tempData = []
pressureData = []
cStarData = []
timeData = []
positionData = []

def graphTime(y, ylabel, yunits):
    #print(f"{ylabel}: {y}\n")
    filename = ylabel + "VsResidenceTime.png"
    plt.figure()
    plt.plot(list(range(timeSteps)), y)
    plt.xlabel('Time (s))')
    plt.ylabel(f"{ylabel} ({yunits})")
    plt.savefig(script_dir / filename)
    plt.close()

def graphPos(y, ylabel, yunits):
    #print(f"{ylabel}: {y}\n")
    filename = ylabel + "VsPosition.png"
    plt.figure()
    plt.plot(list(range(steps)), y)
    plt.xlabel('Position (m)')
    plt.ylabel(f"{ylabel} ({yunits})")
    plt.savefig(script_dir / filename)
    plt.close()
    

def cstarefficiency(): #will update  this function (copy from zero dim code) after loop works
    return 

#outer loop iterates over time
for i in range(timeSteps):
    #inner loop gets all values at every position at one exact time
    for j in range(steps):
        #1: setup --> essentially doing the same exact thing as above
        upstream.phase.TDY = reactor1.phase.TDY
        #network.reinitialize()
        #compute velocity and transform into time
        u[j] = mtotal / area / reactor1.phase.density
        t_reactor1[j] = reactor1.mass / mtotal  # residence time in this reactor
        t[j] = np.sum(t_reactor1)
        #write output data
        states.append(reactor1.phase.state) #at the end of the loop, gives array of the state at every position at this exact time 
    print(states)
    graphPos(states.T, "Temperature", "K")
    mixture = ct.Solution("nDodecane_Reitz.yaml")
    mixture.TPX = reactor1.phase.TPX
    states = ct.SolutionArray(mixture)
    network.advance(network.time + dt)

#current state of this loop: the temp is changing over time, but not across reactors --> probably problem with line 105 or smth in inner loop? 
#ie. no reaction across reactors
#possible fix would be to rewrite the inner loop to resemble the setup before the loop, creating a new reactor each time
#maybe switch structure so time is inner and position is outer? --> would that just cause the same issue?




#graph outputs
#graphTime(tempData, "Temperature", "K")
#graphPos(tempData, "Temperature", "K")

#graphTime(pressureData, "Pressure", "Pa")
#graphPos(pressureData, "Pressure", "Pa")
