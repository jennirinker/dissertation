"""
fixing MO length calculations
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import scipy.io as scio
import numpy as np

dataset = 'NREL'
#fname = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
#    'processed_data\\NREL-metadata.mat'
#fields, py = jr.loadmetadata(fname)
#clean = jr.screenmetadata(fields,py,'NREL')
rec = (2012, 2, 16, 6, 0)
ht = 30

fpath = jr.time2fpath(dataset,rec)
struc20 = scio.loadmat(fpath)

# get all the time series we need
loht_DP, hiht_DP = jr.interpolationHeights(dataset,ht,'Dewpt_Temp')
outdict  = jr.loadtimeseries(dataset,'Dewpt_Temp',loht_DP,struc20)
DPt_lo = outdict['clean']    # celsius
outdict  = jr.loadtimeseries(dataset,'Dewpt_Temp',hiht_DP,struc20)
DPt_hi = outdict['clean']    # celsius
outdict  = jr.loadtimeseries(dataset,'Dewpt_Temp',3,struc20)
DP0t = outdict['clean']    # celsius
outdict = jr.loadtimeseries(dataset,
                        'Temperature',3,struc20) # celsius
T0t = outdict['clean']
outdict = jr.loadtimeseries(dataset,'Pressure',3,struc20)
P0t = outdict['clean']
loht_T, hiht_T = jr.interpolationHeights(dataset,ht,'Temperature')
outdict  = jr.loadtimeseries(dataset,'Temperature',loht_T,struc20)
Tt_lo = outdict['clean']    # celsius
outdict  = jr.loadtimeseries(dataset,'Temperature',hiht_T,struc20)
Tt_hi = outdict['clean']    # celsius

# ================= calculate values at reference height ======================
g, R = 9.81, 287

# get mean values from time histories
DP0 = np.nanmean(DP0t)
P0  = np.nanmean(P0t)
T0  = np.nanmean(T0t)
T0_K = T0 + 273.15

# get derived values
if (DP0 > 0): A, B = 7.5, 237.3
else:            A, B = 9.5, 265.5
e0 = 6.11 * 10 ** ((DP0*A)/(DP0 + B))
q0 = e0 / P0
Tv0 = (T0_K)*(1 + 0.61*q0)
dPdz = - (g * P0) / (R * Tv0)

# ================= calculate values at measurement height ====================

# get mean values from time histories
DP_lo = np.nanmean(DPt_lo)
DP_hi = np.nanmean(DPt_hi)
T_lo  = np.nanmean(Tt_lo)
T_hi  = np.nanmean(Tt_hi)

# interpolate values
DPz = jr.interpolateparameter(dataset,ht,np.nanmean(DP_lo),
                                 np.nanmean(DP_hi),'Dewpt_Temp')
Tz  = jr.interpolateparameter(dataset,ht,np.nanmean(T_lo),
                                 np.nanmean(T_hi),'Temperature')
Tz_K = Tz + 273.15

# get derived values
if (DPz > 0): A, B = 7.5, 237.3
else:            A, B = 9.5, 265.5
ez = 6.11 * 10 ** ((DPz*A)/(DPz + B))
Pz = P0 + (ht - 3)*dPdz
qz = ez / Pz
Tvz_K = (Tz_K)*(1 + 0.61*qz)
Tvz = Tvz_K - 273.15

print(DPz,Pz,Tz,ez,qz)

print('Vapor pressure at height:   {:.2f}'.format(ez))
print('Spec. humidity at height:   {:.3f}'.format(qz))
print('Temperature at height:      {:.2f}'.format(Tz))
print('Virtual temp. at height:    {:.2f}'.format(Tvz))

