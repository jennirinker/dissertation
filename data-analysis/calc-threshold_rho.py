"""
Calculate the threshold values for different sampling frequencies
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
import scipy.stats

Q = 0.99
fs = [10.,20.,50.]
T  = 600.

# -----------------------------------------------------------------------------

for f in fs:
    n_t  = T*f
    n_f  = jr.uniqueComponents(n_t)
    n_dt = n_f - 1
    rhoThresh = np.sqrt(scipy.stats.chi2.ppf(Q,2)/2./n_dt)
    print(n_dt,'{:.3f}'.format(rhoThresh))