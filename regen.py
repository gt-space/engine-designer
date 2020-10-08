from engine_designer.dataCollectionScript import Engine
from engine_designer.regen.regenAnalysis import regenJacket
from engine_designer.regen.regenAnalysis_rat import design_regen_rat
import matplotlib.pyplot as plt

# (thrust (N), chamber pressure (bar), contraciton ratio)
masses = []
coolant_temps = []
combination = []
# a = 5000
for i in range(10):
    try:
        engine = Engine(13000, 80, 2+i) # Create engine object
        engine.design_engine() # Run engine design procedures
        jacket = regenJacket(engine) # Create jacket object
        (profile, T_co, mass) = jacket.get_geometry() # Generate channel geometry
        masses.append(mass)
        coolant_temps.append(T_co)
        # combination.append(a * mass + T_co)
        # if i > 0:
        #     if (coolant_temps[-2] - coolant_temps[-1])/coolant_temps[-1] < 0.03 :
        #         print("CONTRACTION RATIO: " + str(2+i))
        #         break
    except:
        break
plt.plot(coolant_temps)
plt.show()
