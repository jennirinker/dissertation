"""
Demosntrating the differences between removing the slope and matching the end 
points of a time series on the calculation of rho.
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np

# load dataset
dataset          = 'NREL'
fpath            = jr.metadataFpath(dataset)
fields, metadata = jr.loadmetadata(fpath)
clean            = jr.screenmetadata(fields,metadata,dataset)

n_recs = 200

rhos_det = np.empty(n_recs)
rhos_end = np.empty(n_recs)
rhos_not = np.empty(n_recs)
# loop through records, caluclating rho
for i in range(n_recs):
    
    # choose and load random time series
    rec_idx = np.random.choice(clean.shape[0])
    rec     = clean[rec_idx,fields.index('Record_Time')]
    ht      = clean[rec_idx,fields.index('Height')]
    outdict = jr.loadtimeseries(dataset,'Sonic_u',15,rec)
    
    t    = outdict['time']
    u    = outdict['raw']
    u_dt = outdict['clean']
    
    u_ns = jr.remove_spikes(u)[0]
    du   = u_ns[-1] - u_ns[0]
    u_nd = u_ns - du/(t.size*t[1])*t
    u_nd = u_nd - np.mean(u_nd) + np.mean(u_ns)
    
    rhos_det[i] = jr.signalPhaseCoherence(u_dt)[0]
    rhos_end[i] = jr.signalPhaseCoherence(u_nd)[0]
    rhos_not[i] = jr.signalPhaseCoherence(u_ns)[0]
    
    
#print('Detrend: {:.3f}'.format(rho_det))
#print('End pts: {:.3f}'.format(rho_end))
#print('Diff:    {:.3f}'.format(rho_det-rho_end))
#
#plt.figure(1)
#plt.clf()
#plt.plot(t,u)
#plt.plot(t,u_dt)
#plt.plot(t,u_nd)
