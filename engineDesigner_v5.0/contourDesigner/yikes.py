import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

all_species = {S.name: S for S in ct.Species.list_from_file('CombProducts.yaml',section='species')}
JetA_X = {
    "C10H22(N)":  0.85,   # n-decane
    "C3H8":     0.15    # propane
}

LOX = ct.Solution(thermo='ideal-gas',species=[all_species['O2']])
JetA = ct.Solution(thermo='ideal-gas',species=[all_species[S] for S in JetA_X])

Ts=np.linspace(2000,4000,2000)
for name in all_species.keys():
    spec=all_species[name]
    enth=[]
    for T in Ts:
        sol = ct.Solution(thermo='ideal-gas',species=[spec])
        sol.TP=T,1e5
        sol.equilibrate('TP')
        enth.append(sol.entropy_mass)

    fig, ax= plt.subplots()
    ax.plot(Ts, enth)
    ax.set_title(name)
    plt.show()

print(len(all_species.keys()))