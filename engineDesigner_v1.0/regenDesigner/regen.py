import sys
import math
import numpy as np
import matplotlib.pyplot as plt
# from .bartzCorrelation import Bartz
# from .jetAProps import getProps
# from .pressureDrop import pressureDrop

np.set_printoptions(threshold=sys.maxsize) # Print full arrays (for debugging)

class regenJacket:
    # Initialize the jacket object upon declaration
    def __init__(self, engine, channel_h=0.001, wall_t=0.001, min_fin_w = 0.001, min_channel_w = 0.0015875/2, T_wg=750, T_co=462):
        self.engine = engine # Engine object to be jacketed
        self.channel_h = channel_h # Channel height (m)
        self.wall_t = wall_t # Inner wall thickness (m)
        self.min_fin_w = min_fin_w # Minimum manufacutrable fin width (m)
        self.min_channel_w = min_channel_w # Minimum manufacutrable channel width (m)
        self.T_wg = T_wg # Gas-side wall temperature initial estimate (K)
        self.T_co = T_co # Desired coolant outlet temperature (K)
        self.T_max = 850 # Max allowed temp for material (copper melts at 1358 K)
