from engine_designer.dataCollectionScript import Engine
from engine_designer.regen.regenAnalysis import regenJacket
from engine_designer.regen.regenAnalysis_rat import design_regen_rat
from conrat_solver import get_conrat
import matplotlib.pyplot as plt
import numpy as np

# Inputs
thrust = 17792.89 # thrust (N)
C_p = 48.2633 # chamber pressure (bar)
conrat = 6
T_co_max = 462 # Max allowed coolant outlet temp (K)
L_star_max = 1.27 # Max allowed L-star (m) (H&H pg. 72)

def get_conrat(thrust, P_c, plot = False):
    # Solve for contraction ratio that first meets desired coolant outlet temp
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
            # This may error as very high contraction ratios are not geometrically solvable
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

def l_star():
    L_star= 1 # Lowest acceptable characteristic length (m)
    for i in range(12):
        L_star += i * 0.01 # Increment L-star
        engine = Engine(thrust, C_p, conrat, LStar = L_star) # Create engine object
        engine.design_engine() # Run engine design procedures
        jacket = regenJacket(engine) # Create jacket object
        (profile, T_co, mass) = jacket.get_geometry() # Generate channel geometry
        print(T_co)
    return L_star

def regen(L_star):
    # conrat = get_conrat(thrust, C_p, True) # Run this to solve for conrat based on minimizing outlet temp (takes a while)
    engine = Engine(thrust, C_p, conrat, LStar = L_star) # Create engine object
    engine.design_engine() # Run engine design procedures
    jacket = regenJacket(engine) # Create jacket object
    (profile, T_co, mass) = jacket.get_geometry() # Generate channel geometry

if __name__ == "__main__":
    L_star = l_star()
    regen(L_star)
