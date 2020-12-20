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
MR = 2.0 # Mixture ratio (lox to fuel)
T_co_max = 462 # Max allowed coolant outlet temp (K)
L_star_max = 1.27 # Max allowed L-star (m) (H&H pg. 72)

def get_modeling_data(engine, profile, num_channels):
    # Display relevant data for modeling
    print(" ---=== DATA FOR MODELING ===--- ")
    print("All dimensions in meters")
    print("Total Length: ", engine.engineProps[-1, 1])
    print("Injector Radius: ", engine.engineProps[0, 0])
    print("Throat Radius: ", engine.R_t)
    print("Exit Radius: ", engine.engineProps[-1, 0])
    print("Throat Contour Radius: ", engine.R_t * 1.5)
    print("Lead in Radius: ", engine.conLeadInRadius)
    print("Number of Channels: ", num_channels)
    print("Channel Barrel Width: ", profile[0, 1])
    print("Channel Throat Width: ", profile[engine.throatInd + 100, 1])
    print("Channel Exit Width: ", profile[-1, 1])

def get_conrat(thrust, P_c, plot = False):
    # Solve for contraction ratio that first meets desired coolant outlet temp
    coolant_temps = []
    outlet = 400 # Max allowed coolant outlet temperature
    conrats = []
    masses = []

    for i in range(20):
        try:
            conrat = 2+i/2
            engine = Engine(thrust, P_c, conrat, MR = MR) # Create engine object
            engine.design_engine() # Run engine design procedures
            jacket = regenJacket(engine, T_co=outlet) # Create jacket object
            (profile, T_co, mass, num_channels) = jacket.get_geometry() # Generate channel geometry
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
        engine = Engine(thrust, C_p, conrat, LStar = L_star, MR = MR) # Create engine object
        engine.design_engine() # Run engine design procedures
        print(engine.engineProps)
        jacket = regenJacket(engine) # Create jacket object
        (profile, T_co, mass, num_channels) = jacket.get_geometry() # Generate channel geometry
        if (T_co > T_co_max) or L_star >= L_star_max:
            break
    return L_star

def regen(L_star):
    # conrat = get_conrat(thrust, C_p, True) # Run this to solve for conrat based on minimizing outlet temp (takes a while)
    engine = Engine(thrust, C_p, conrat, LStar = L_star, MR = MR) # Create engine object
    engine.design_engine() # Run engine design procedures
    jacket = regenJacket(engine) # Create jacket object
    (profile, T_co, mass, num_channels) = jacket.get_geometry() # Generate channel geometry
    print(profile)
    get_modeling_data(engine, profile, num_channels)

if __name__ == "__main__":
    L_star = l_star()
    regen(L_star)
