"""
Debugging fatigue DEL calculations
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\gist\\rainflow'
if (libpath not in sys.path): sys.path.append(libpath)

import numpy as np
import matplotlib.pyplot as plt
from WISDEM_rainflow.rainflow import determine_peaks
import JR_Library.main as jr
from JR_Library.peakdetect import peakdetect
from rainflow import rainflow

# choose whether to reload FAST dictionary
reload_dict = 0

# simulation parameters
t_lookahead = 0.1           # time to look ahead for gist algorithm [s]
key         = 'RootMOoP1'   # what value to analyze
m           = 10.           # SN slope
equivFreq   = 1.0           # equivalent frequency [Hz]
T_DEL       = 600.          # time span for DEL calculation [s]
PSF         = 1.            # partial safety factor

# path to FAST file
fpath_fast = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation\\' + \
                'FAST_models\\FAST7\\WP0.75A08V00_equil\\WP0.75A08V00_24134.out'

# reload dictionary if requested
if reload_dict:
    fast_dict = jr.ReadFASTFile(fpath_fast)
    
# get time series from dictionary
data = fast_dict['Data']
fields = fast_dict['Fields']
units = fast_dict['Units']
t         = data[:,fields.index('Time')]
y         = data[:,fields.index(key)]
key_units = units[fields.index(key)]

# -------------------------- get extreme values -------------------------------

# number of steps to look ahead for Gist algorithm
lookahead = int(t_lookahead / t[1])

# get peak values and indices using WISDEM algorithm
peaks, ipeaks = determine_peaks(y)
xpeaks_W = t[ipeaks]
ypeaks_W = y[ipeaks]

# get peak values and indices using gist algorithm
maxinfo, mininfo = peakdetect(y, t, 
                              lookahead=1)
xpeaks_g = np.array([tup[0] for tup in maxinfo] + [tup[0] for tup in mininfo])
ypeaks_g = np.array([tup[1] for tup in maxinfo] + [tup[1] for tup in mininfo])

# resort gist algorithm into increasing x value
isort    = xpeaks_g.argsort()
xpeaks_g = xpeaks_g[isort]
ypeaks_g = ypeaks_g[isort]

# -------------------------- rainflow counting --------------------------------

N_eq = T_DEL * equivFreq

rainflow_W = rainflow(ypeaks_W)
rainflow_g = rainflow(ypeaks_g)

n_W, S_W = rainflow_W[3,:], rainflow_W[2,:]
n_g, S_g = rainflow_g[3,:], rainflow_g[2,:]

DEL_W = ( np.sum(n_W * ( S_W )**m) / ( N_eq ) ) ** (1. / m) * PSF
DEL_g = ( np.sum(n_g * ( S_g )**m) / ( N_eq ) ) ** (1. / m) * PSF
perc_diff = ((DEL_g - DEL_W)/DEL_W)*100.

print('\n10-minute DELs for {:s}:'.format(key))
print('--------------------------------')
print('WISDEM peak detection: {:8.1f} {:s}'.format(DEL_W,key_units))
print('Gist peak detection:   {:8.1f} {:s}'.format(DEL_g,key_units))
print('\nPercent difference:   {:8.1f}%'.format(perc_diff))



