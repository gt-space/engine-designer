from engine_designer.dataCollectionScript import Engine
from engine_designer.regen.regenAnalysis import regenJacket
from engine_designer.regen.regenAnalysis_rat import design_regen_rat
from conrat_solver import get_conrat
import matplotlib.pyplot as plt

thrust = 13000 # thrust (N)
C_p = 30 # chamber pressure (bar)
conrat = 6
# conrat = get_conrat(thrust, C_p) # Run this to solve for conrat based on minimizing outlet temp (takes a while)
engine = Engine(thrust, C_p, conrat) # Create engine object
engine.design_engine() # Run engine design procedures
jacket = regenJacket(engine) # Create jacket object
(profile, T_co, mass) = jacket.get_geometry() # Generate channel geometry
