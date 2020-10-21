from engine_designer.dataCollectionScript import Engine
from engine_designer.regen.regenAnalysis import regenJacket
import matplotlib.pyplot as plt
import numpy as np

# Solve for contraction ratio that first meets desired coolant outlet temp

def get_conrat(thrust, P_c, plot = False):
    # (thrust (N), chamber pressure (bar), contraciton ratio)
    coolant_temps = []
    outlet = 400 # Max allowed coolant outlet temperature
    conrats = []
    masses = []

    for i in range(20):
        try:
            conrat = 2+i/2
            engine = Engine(thrust, P_c, conrat) # Create engine object
            engine.design_engine() # Run engine design procedures
            jacket = regenJacket(engine, T_co=outlet) # Create jacket object
            (profile, T_co, mass) = jacket.get_geometry() # Generate channel geometry
            coolant_temps.append(T_co)
            conrats.append(conrat)
            masses.append(mass)
            if coolant_temps[i] < outlet:
                break
        except:
            break
    if plot:
        print('CONTRACTION RATIO: ' + str(conrat))
        plt.subplot(1, 2, 1)
        plt.plot(conrats, coolant_temps)
        plt.ylabel('Coolant Outlet Temp (K)', fontsize=16)
        plt.xlabel('Contraction Ratio', fontsize=16)

        plt.subplot(1, 2, 2)
        plt.plot(conrats, masses)
        plt.ylabel('Mass (kg)', fontsize=16)
        plt.xlabel('Contraction Ratio', fontsize=16)
        plt.show()

    return conrat
