"""
Messing around with writing .bat files for Monsoon calculations
"""

import numpy as np

# simulation input variables
IndSim  = 1
NSims   = 6000


# dervied variables
NDigits = int(np.ceil(np.log10(NSims)))
SimID   = '{:d}'.format(IndSim).zfill(NDigits)
SimName = 'Bat{:s}'.format(SimID)

