"""
Debug NREL processing routine in Python
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
import scipy.io as scio
import scipy.signal
import calendar, time

# %% ===================== load metadata structures ===========================

## load fields, metadata from matlab
#flds_mat, raw_mat = jr.loadNRELmatlab()
#clean_mat = jr.screenmetadata(flds_mat,raw_mat,'NREL')
#
## load python-processed metadata
fname = 'C:\\Users\\jrinker\\Dropbox\\research' + \
    '\\processed_data\\NREL-metadata.mat'
flds_py,raw_py = jr.loadmetadata(fname)
clean_py =  jr.screenmetadata(flds_py,raw_py,'NREL')
#del raw_mat, raw_py

# %% ====================== plot highest rho value ===========================

# get python and matlab indices
#idx_py = 1                                      # python index
idx_py = clean_py[:,
                  flds_py.index('Concentration_u')].argmax()
time_flt = clean_py[idx_py,0]                   # float of timestamp
time_tup = jr.timeflt2tup(time_flt)             # tuple of timestamp
ht       = int(clean_py[idx_py,2])             # measurment height
rec_vec  = np.asarray(time_tup + (ht,))         # time/height tuple

#outdict = jr.loadtimeseries('NREL','Sonic_u',ht,time_tup)
#
#t = outdict['time']
#u_raw = outdict['raw']
#u_cl = outdict['clean']
#
#plt.figure(1,figsize=(6,2.5))
#plt.clf()
#
#plt.plot(t,u_raw)
#plt.plot(t,u_cl)
#plt.xlabel('Time [s]')
#plt.ylabel('Sonic u [m/s]')

