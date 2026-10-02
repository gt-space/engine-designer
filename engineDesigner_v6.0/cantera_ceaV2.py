import cantera as ct
import math
import matplotlib.pyplot as plt
import numpy as np
import time
from pathlib import Path
from rocketcea.cea_obj import CEA_Obj

# initializing cea for rp-1 lox
cea_obj = CEA_Obj(oxName='LOX', fuelName='RP1')

#in current folder, creates a folder for the graphs. graphs vs. time will be placed directly into this folder.
#creates a subfolder in the graphs folder for the position. this folder has two subfolders, one for temperature and one for pressure
#since a new temperature vs. position and a new pressure vs. position graph is created at certain times, the total amount of graphs will be
#equal to the number of timesteps divided by 0.05, doubled plus the 8 vs. time graphs.
script_dir = Path(__file__).parent
graph_dir = script_dir / "cantera_cea_graphs"
graph_dir.mkdir(exist_ok = True)

#conversions
psiaToPa = 6894.76
inch_to_meter = 0.0254

#variables
simTime = 0.15  # 300 ms 
#dt = 5e-4 don't need dt anymore since we are using adaptive time stepping

#initial conditions
mdotFuel = 1.55 #from regen qual trr
mdotLOX = 3.1 #from regen qual trr
mdotTotal = mdotFuel + mdotLOX
initialPc = 250 * psiaToPa #from teams messages
of_ratio = mdotLOX / mdotFuel

#pull official NASA CEA benchmark data
cea_tc = cea_obj.get_Tcomb(Pc=250.0, MR=of_ratio) / 1.8 #convert to kelvin
cea_cstar = cea_obj.get_Cstar(Pc=250.0, MR=of_ratio) * 0.3048  #convert to m/s
print(f"\n[CEA Benchmark] T_comb: {cea_tc:.1f} K | C*: {cea_cstar:.1f} m/s")

#calculate volume of each reactor
x = [] #used to get liner geometry
y = [] #used to get liner geometry
with open(script_dir / "liner_contour_inner_CR5_PE8_MR2.txt", "r") as f: #reads the contour text file where contour is in meters
    for line in f:
        y.append(float(line[0:8]))
        x.append(float(line[9:17]))
#changed data parsing to pull exact contour values instead of polyfit curves
throatRadius = min(y) #throat radius is taken to be smallest radius in the contour
throatArea = math.pi * throatRadius * throatRadius #calculates throat area for cstar calculation
throatIdx = y.index(throatRadius)
chamberRadius = max(y[:throatIdx+1])
convTol = 1e-4
startConv = max(i for i in range(0, throatIdx + 1) if y[i] >= chamberRadius * (1 - convTol))
x_c = np.array(x[:startConv+1])
y_c = np.array(y[:startConv+1])
chamberVolume = float(np.sum(0.5 * np.pi * (y_c[1:]**2 + y_c[:-1]**2) * np.diff(x_c)))

print(f"Volume: {chamberVolume*1000:.2f} L | Throat Area: {throatArea*10000:.2f} cm^2")

# throat area required for target chamber pressure since its all coupled
cstar_for_sizing = cea_cstar 
throatArea_required = mdotTotal * cstar_for_sizing / initialPc
Pc_with_contour = mdotTotal * cstar_for_sizing / throatArea / psiaToPa


#can switch to true to use the actual contour throat area instead
USE_CONTOUR_THROAT = False
#force throat area to match the boundary conditions in case contour doesn't match
#not entirely sure if this is actually needed but a forum post on cantera used the same approach so I just copied it
if not USE_CONTOUR_THROAT:
    throatArea = throatArea_required

#reactor network initialization
#chaining well stirred reactors to try and simulate a plug flow reactor
#artificial volume scaler to change residence time
volume_scaler = 1.0
numReactors = 6 

seg_vol = 0.5 * np.pi * (y_c[1:]**2 + y_c[:-1]**2) * np.diff(x_c)
cum_vol = np.concatenate(([0.0], np.cumsum(seg_vol)))
edges = np.linspace(x_c[0], x_c[-1], numReactors + 1)      # equal-length reactors along x
reactorVolumes = np.diff(np.interp(edges, x_c, cum_vol)) * volume_scaler

reactorList = [] #organizational list of all reactors for easy access

#initialize the gas state for the reactors
gas_init = ct.Solution("nDodecane_Reitz.yaml", "nDodecane_RK")
gas_init.TP = 3000.0, ct.one_atm 
gas_init.set_mixture_fraction(mdotFuel / mdotTotal, fuel="c12h26:1.0", oxidizer="o2:1.0")
gas_init.equilibrate("TP")

#initialize chamber reactors
for i in range(0, numReactors):
    currentReactor = ct.Reactor(gas_init, volume=reactorVolumes[i], clone=True)
    reactorList.append(currentReactor)

#properties to simulate the injectors for more accurate transient start-up behavior
fuelGas = ct.Solution("nDodecane_Reitz.yaml", "nDodecane_RK")
fuelGas.TPX = 298.15, initialPc * 1.5, "c12h26:1.0"
oxGas = ct.Solution("nDodecane_Reitz.yaml", "nDodecane_RK")
oxGas.TPX = 298.15, initialPc * 1.5, "o2:1.0"

#the point of the reservoirs is to provide a constant pressure source for the mass flow controllers to pull from
#this is also from a forum post so I'm not entirely sure if it's actually applicable here
fuelRes = ct.Reservoir(fuelGas, clone=True)
oxRes = ct.Reservoir(oxGas, clone=True)

#trying to simulate valve opening and ramping over 30 ms for a more realistic transient start-up
#completley arbitrary but worked well for this case
t_ramp_duration = 0.03 
def ramp(t):
    return min(t / t_ramp_duration, 1.0)

#temp lambda functions to allow for time-dependent mass flow rates
mfcFuel = ct.MassFlowController(fuelRes, reactorList[0])
mfcFuel.mass_flow_rate = lambda t: mdotFuel * ramp(t)

#mass flow controller for oxidizer
mfcOx = ct.MassFlowController(oxRes, reactorList[0])
mfcOx.mass_flow_rate = lambda t: mdotLOX * ramp(t)

#forcing a pressure drop across the reactors to try and simulate an actual pressure profile
#again, completely arbitrary and was found by just varying the pressure drop until the last reactor was at the target chamber pressure
target_dp_total = 0.02 * initialPc
dp_per_valve = target_dp_total / max(1, numReactors - 1)
derived_Kv = mdotTotal / dp_per_valve

#connect reactors
valves = [] #used to store valve objects
for i in range(1, len(reactorList), 1):
    previousReactor = reactorList[i-1]
    currentReactor = reactorList[i]
    connection = ct.Valve(previousReactor, currentReactor, K=derived_Kv) #each reactor is connected by a valve object, valve constant chosen arbitrarily
    valves.append(connection)

#set up outlet
outletSolution = ct.Solution("nDodecane_Reitz.yaml", "nDodecane_RK")
outletSolution.TPX = 298.15, ct.one_atm, "n2:1.0"
outletRes = ct.Reservoir(outletSolution, clone=True)

def choked_nozzle_mdot(t):
    p = reactorList[-1].phase.P
    #cubic Hermite smoothing to prevent ode integrator noise near the nozzle opening
    x = min(max((p - ct.one_atm) / (0.5 * ct.one_atm), 0.0), 1.0)
    opening = x * x * (3.0 - 2.0 * x)
    if opening == 0.0:
        return 0.0
    T = reactorList[-1].phase.T
    gam = reactorList[-1].phase.cp_mass / reactorList[-1].phase.cv_mass
    R = ct.gas_constant / reactorList[-1].phase.mean_molecular_weight
    psi = math.sqrt(gam) * ((2 / (gam + 1)) ** ((gam + 1) / (2 * (gam - 1))))
    return opening * throatArea * (p / math.sqrt(R * T)) * psi

#mass flow controller to simulate choked nozzle flow at the exit of the last reactor
nozzle_mfc = ct.MassFlowController(reactorList[-1], outletRes, mdot=choked_nozzle_mdot)

#set up reactor network
network = ct.ReactorNet(reactorList)
network.max_err_test_fails = 12  # default is 7, give CVODE more attempts
network.atol = 1e-9
network.rtol = 1e-6

#calculate ideal cstar by creating equilibrium mixture
gas_eq = ct.Solution("nDodecane_Reitz.yaml", "nDodecane_RK")
gas_eq.TP = 298.15, initialPc
gas_eq.set_mixture_fraction(mdotFuel/mdotTotal, fuel="c12h26:1.0", oxidizer="o2:1.0") #mixture ratio comes from vespula requirements
gas_eq.equilibrate("HP")
gamma = gas_eq.cp_mass / gas_eq.cv_mass
R_u = ct.gas_constant
M = gas_eq.mean_molecular_weight
R_specific = R_u / M
cstar_ideal = math.sqrt((1 / gamma) * (((gamma + 1) / 2) ** ((gamma + 1) / (gamma - 1))) * (R_specific * gas_eq.T))

# cea benchmark data for enthalpy deficits of liquid propellants
plot_species = ["co", "co2", "h2", "h2o", "oh", "h", "o", "o2"]
cea_species_key = {"co": "*CO", "co2": "*CO2", "h2": "*H2", "h2o": "H2O",
                   "oh": "*OH", "h": "*H", "o": "*O", "o2": "*O2"}
cea_X, cea_MW, cea_cp = {}, None, None

_, _mf = cea_obj.get_SpeciesMoleFractions(Pc=250.0, MR=of_ratio, eps=4.0)
for sp, key in cea_species_key.items():
    cea_X[sp] = float(_mf[key][1]) if key in _mf else np.nan   # index 1 = chamber
cea_MW = float(cea_obj.get_Chamber_MolWt_gamma(Pc=250.0, MR=of_ratio)[0])

#RocketCEA gives Cp in BTU/(lbm*R) so we convert
#   also includes the heat absorbed by shifting dissociation reactions, not just heating the gas

cea_cp = float(cea_obj.get_Chamber_Cp(Pc=250.0, MR=of_ratio)) * 4186.8
print(f"[RocketCEA] Chamber MW = {cea_MW:.3f} | Cp(equilibrium) = {cea_cp:.1f} J/kg/K")

#have to manually calculate equilibrium Cp since Cantera defaults to frozen Cp and doesn't have a built-in function
#this numerical approximation is called central finite differences

'''
we're basically using a slightly higher and lower temperature to then solving
for a slightly higher and lower enthalpy, then using the difference in enthalpy
over the difference in temperature to get the Cp
'''

gas_scratch = ct.Solution("nDodecane_Reitz.yaml", "nDodecane_RK")
def equilibrium_cp(T, P, X, dT=5.0):
    gas_scratch.TPX = T + dT, P, X
    gas_scratch.equilibrate("TP")
    h_plus = gas_scratch.enthalpy_mass
    gas_scratch.TPX = T - dT, P, X
    gas_scratch.equilibrate("TP")
    h_minus = gas_scratch.enthalpy_mass
    return (h_plus - h_minus) / (2.0 * dT)

#transient simulation
print("\nRunning Simulation...")
start_time = time.time()

#one label per reactor, reused in every legend/title so basically reactor 1 = injector end, last = nozzle end
reactor_labels = [f"Reactor {i+1}" for i in range(numReactors)]
reactor_labels[0] += " (injector)"
reactor_labels[-1] += " (nozzle)"
reactor_colors = plt.cm.viridis(np.linspace(0.0, 0.9, numReactors))

#shared time axis and commanded propellant flow, plus one list per reactor for each recorded quantity
timeArr, mdot_in_arr = [], []
p_r  = [[] for _ in range(numReactors)]     
t_r  = [[] for _ in range(numReactors)]    
mw_r = [[] for _ in range(numReactors)]     
cp_r = [[] for _ in range(numReactors)]     
species_r = [{sp: [] for sp in plot_species} for _ in range(numReactors)]   #mole fractions

current_time = 0.0
step_count = 0

#adaptive time stepping loop to handle the transient simulation with varying time steps based on the current time
while current_time < simTime:
    if current_time < 0.01:
        dt_step = 1e-5
    elif current_time < 0.05:
        dt_step = 1e-4
    else:
        dt_step = 1e-3

    current_time = min(simTime, current_time + dt_step)

    try:
        network.advance(current_time)
    except Exception as e:
        print(f"\nSolver failed at t = {network.time:.5f}s: {e}")
        break

    #record data for each reactor
    t = network.time
    timeArr.append(t) #adds the current timestamp to the array
    mdot_in_arr.append(mdotTotal * ramp(t))

    #equilibrium Cp tanks runtime so we only calculate it once every 6 steps after the first 20 steps
    do_cp = (t > 0.02 and step_count % 6 == 0)

    for i, r in enumerate(reactorList):
        phi = r.phase                       # live state of reactor i
        p_r[i].append(phi.P / psiaToPa)
        t_r[i].append(phi.T)
        mw_r[i].append(phi.mean_molecular_weight)
        for sp in plot_species:
            species_r[i][sp].append(phi.X[phi.species_index(sp)])
        if do_cp:
            try:
                cp_r[i].append(equilibrium_cp(phi.T, phi.P, phi.X))
            except Exception:
                cp_r[i].append(np.nan)
        else:
            cp_r[i].append(np.nan)

    step_count += 1
    if step_count % 50 == 0:
        print(f"t = {t:.4f}s | P_R1 = {p_r[0][-1]:6.1f} psia | P_R{numReactors} = {p_r[-1][-1]:6.1f} psia | T_R1 = {t_r[0][-1]:6.1f} K | T_R{numReactors} = {t_r[-1][-1]:6.1f} K")

print(f"\nSimulation Complete in {time.time() - start_time:.2f} seconds.")

#post processing
time_np = np.asarray(timeArr)
mdot_in_np = np.asarray(mdot_in_arr)
p_np  = [np.asarray(x) for x in p_r]
t_np  = [np.asarray(x) for x in t_r]
mw_np = [np.asarray(x) for x in mw_r]
cp_np = [np.asarray(x) for x in cp_r]

# measured C-star = P * At / mdot_in, evaluated with each reactor's own pressure
valid = mdot_in_np >= 0.5 * mdotTotal

# evaluate C-star at a range of pressures using RocketCEA for comparison
p_all = np.concatenate(p_np)
p_grid = np.linspace(10.0, 1.05 * p_all.max(), 40)
c_grid = np.array([cea_obj.get_Cstar(Pc=p, MR=of_ratio) * 0.3048 for p in p_grid])

cstar_meas_r, cstar_ref_r, cstar_eff_r = [], [], []
for i in range(numReactors):
    meas = np.where(valid, p_np[i] * psiaToPa * throatArea / np.maximum(mdot_in_np, 1e-12), np.nan)
    ref = np.interp(p_np[i], p_grid, c_grid)
    cstar_meas_r.append(meas)
    cstar_ref_r.append(ref)
    cstar_eff_r.append(100.0 * meas / ref)

#plots
# 1. Chamber Composition vs CEA: one panel per reactor
ncols = 3
nrows = math.ceil(numReactors / ncols)
fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), sharex=True, sharey=True)
axes = np.atleast_1d(axes).ravel()
sp_colors = plt.cm.tab10(np.linspace(0, 1, len(plot_species)))
for i in range(numReactors):
    ax_i = axes[i]
    for c, sp in zip(sp_colors, plot_species):
        ax_i.plot(time_np, species_r[i][sp], color=c, linewidth=1.8, label=sp.upper())
        if np.isfinite(cea_X.get(sp, np.nan)):
            ax_i.axhline(cea_X[sp], color=c, linestyle='--', linewidth=1.0)
    ax_i.set_yscale("log")
    ax_i.set_ylim(1e-4, 1.0)
    ax_i.set_title(reactor_labels[i])
    ax_i.grid(True, which="both", linestyle='--', alpha=0.4)
    if i % ncols == 0:
        ax_i.set_ylabel("Mole Fraction [-]")
    if i >= numReactors - ncols:
        ax_i.set_xlabel("Time [s]")
for j in range(numReactors, len(axes)):      # hide unused panels if reactors don't fill the grid
    axes[j].set_visible(False)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=len(plot_species), fontsize=9)
fig.suptitle("Composition by Reactor (solid = Cantera, dashed = RocketCEA equilibrium)")
plt.tight_layout(rect=[0, 0.05, 1, 0.96])
plt.savefig(graph_dir / "1_Chamber_Composition_vs_CEA.png", dpi=150)
plt.close()

# 2. Molecular Weight and Cp vs CEA: one line per reactor
fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
for i in range(numReactors):
    ax[0].plot(time_np, mw_np[i], color=reactor_colors[i], linewidth=1.8, label=reactor_labels[i])
ax[0].axhline(cea_MW, color='purple', linestyle='--', label=f"RocketCEA ({cea_MW:.2f})")
ax[0].axhline(M, color='black', linestyle=':', label=f"Cantera equilibrium ({M:.2f})")
ax[0].set_ylabel("Mean Molecular Weight [g/mol]")
ax[0].set_title("Chamber MW and Cp vs RocketCEA (by reactor)")
ax[0].grid(True, linestyle='--', alpha=0.6)
ax[0].legend(fontsize=8, ncol=2)

# Equilibrium Cp only (comparable to CEA); frozen cp_mass is ~half as large and not comparable
cp_all = []
for i in range(numReactors):
    okc = np.isfinite(cp_np[i])
    ax[1].plot(time_np[okc], cp_np[i][okc], 'o-', color=reactor_colors[i], markersize=3,
               linewidth=1.5, label=reactor_labels[i])
    cp_all.extend(cp_np[i][okc])
ax[1].axhline(cea_cp, color='purple', linestyle='--', label=f"RocketCEA ({cea_cp:.0f} J/kg/K)")
if len(cp_all) > 0:
    y_min = min(np.min(cp_all), cea_cp)
    y_max = max(np.max(cp_all), cea_cp)
    pad = max(50.0, 0.1 * (y_max - y_min))
    ax[1].set_ylim(y_min - pad, y_max + pad)
ax[1].set_xlabel("Time [s]")
ax[1].set_ylabel(r"Equilibrium $c_p$ [J/kg/K]")
ax[1].grid(True, linestyle='--', alpha=0.6)
ax[1].legend(fontsize=8, ncol=2)
plt.tight_layout()
plt.savefig(graph_dir / "2_MW_Cp_vs_CEA.png", dpi=150)
plt.close()

# 3. Temperature by reactor vs CEA
plt.figure(figsize=(9, 5))
for i in range(numReactors):
    plt.plot(time_np, t_np[i], color=reactor_colors[i], linewidth=1.8, label=reactor_labels[i])
plt.axhline(gas_eq.T, color='black', linestyle=':', label=f"Cantera Equilibrium ({gas_eq.T:.0f} K)")
plt.axhline(cea_tc, color='purple', linestyle='--', label=f"RocketCEA T_comb ({cea_tc:.0f} K)")
plt.xlabel("Time [s]")
plt.ylabel("Temperature [K]")
plt.title("Temperature by Reactor vs Equilibrium Benchmarks")
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=8, ncol=2)
plt.tight_layout()
plt.savefig(graph_dir / "3_Chamber_Temperature_vs_CEA.png", dpi=150)
plt.close()

# 4. C* and Efficiency by reactor vs CEA
fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
for i in range(numReactors):
    ax[0].plot(time_np, cstar_meas_r[i], color=reactor_colors[i], linewidth=1.8, label=reactor_labels[i])
ax[0].plot(time_np, cstar_ref_r[-1], color='purple', linestyle='--', label="RocketCEA C* (at nozzle-end Pc)")
ax[0].axhline(cstar_ideal, color='black', linestyle=':', label=f"Cantera equilibrium ({cstar_ideal:.0f} m/s)")
ax[0].set_ylabel("C* [m/s]")
ax[0].set_title(r"Characteristic Velocity vs RocketCEA ($C^*=P_iA_t/\dot m_{in}$ per reactor)")
ax[0].grid(True, linestyle='--', alpha=0.6)
ax[0].legend(loc="lower right", fontsize=8, ncol=2)

eff_all = []
for i in range(numReactors):
    ax[1].plot(time_np, cstar_eff_r[i], color=reactor_colors[i], linewidth=1.8, label=reactor_labels[i])
    eff_all.extend(cstar_eff_r[i][np.isfinite(cstar_eff_r[i])])
ax[1].axhline(100.0, color='purple', linestyle='--')
if len(eff_all) > 0:
    lo, hi = np.min(eff_all), np.max(eff_all)
    pad = max(1.0, 0.1 * (hi - lo))
    ax[1].set_ylim(lo - pad, hi + pad)
ax[1].set_xlabel("Time [s]")
ax[1].set_ylabel(r"$\eta_{C^*}$ [%]")
ax[1].grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig(graph_dir / "4_Cstar_vs_CEA.png", dpi=150)
plt.close()

#summary table (nozzle-end reactor = the true chamber/stagnation state, compared to CEA)
last = slice(-5, None)
cp_eq_final = np.nanmean(cp_np[-1][-30:]) if np.any(np.isfinite(cp_np[-1])) else np.nan
print(f"\n=========== {reactor_labels[-1]} vs rocketcea ===========")
print(f"{'Quantity':<26}{'Cantera':>14}{'RocketCEA':>14}")
rows = [("Pc [psia]", np.mean(p_np[-1][last]), 250.0),
        ("T_chamber [K]", np.mean(t_np[-1][last]), cea_tc),
        ("MW [g/mol]", np.mean(mw_np[-1][last]), cea_MW),
        ("Cp equilibrium [J/kg/K]", cp_eq_final, cea_cp),
        ("C* [m/s]", np.nanmean(cstar_meas_r[-1][last]), cea_cstar)]
for sp in plot_species:
    rows.append((f"X_{sp.upper()}", np.mean(species_r[-1][sp][-5:]), cea_X.get(sp, np.nan)))
for name, a, b in rows:
    print(f"{name:<26}{a:>14.4f}{b:>14.4f}")

print("\n=========== final state by reactor ===========")
print(f"{'Reactor':<24}{'P [psia]':>10}{'T [K]':>10}{'MW':>9}{'Cp_eq':>10}{'C* [m/s]':>11}")
for i in range(numReactors):
    cp_i = np.nanmean(cp_np[i][-30:]) if np.any(np.isfinite(cp_np[i])) else np.nan
    print(f"{reactor_labels[i]:<24}{np.mean(p_np[i][last]):>10.1f}{np.mean(t_np[i][last]):>10.1f}"
          f"{np.mean(mw_np[i][last]):>9.2f}{cp_i:>10.1f}{np.nanmean(cstar_meas_r[i][last]):>11.1f}")
print(f"\nGraphs saved in: {graph_dir}")