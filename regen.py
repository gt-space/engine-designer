from engine_designer.dataCollectionScript import Engine
from engine_designer.regen.regenAnalysis import regenJacket
from engine_designer.regen.regenAnalysis_rat import design_regen_rat

# (thrust (N), chamber pressure (bar), contraciton ratio)
engine = Engine(3500, 18, 6) # Create engine object
engine.design_engine()
jacket = regenJacket(engine)
jacket.get_geometry()
