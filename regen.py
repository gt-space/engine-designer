from engine_designer.dataCollectionScript import Engine
from engine_designer.regen.regenAnalysis import regenJacket
from engine_designer.regen.regenAnalysis_rat import design_regen_rat
import matplotlib.pyplot as plt

# (thrust (N), chamber pressure (bar), contraciton ratio)
engine = Engine(13000, 30, 6) # Create engine object
engine.design_engine() # Run engine design procedures
jacket = regenJacket(engine) # Create jacket object
(profile, T_co, mass) = jacket.get_geometry() # Generate channel geometry
