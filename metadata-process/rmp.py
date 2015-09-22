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
DP_lo = outdict['clean']    # celsius
outdict  = jr.loadtimeseries(dataset,'Dewpt_Temp',hiht_DP,struc20)
DP_hi = outdict['clean']    # celsius
outdict  = jr.loadtimeseries(dataset,'Dewpt_Temp',3,struc20)
DP0 = outdict['clean']    # celsius
outdict = jr.loadtimeseries(dataset,
                        'Temperature',3,struc20) # celsius
T0 = outdict['clean']
outdict = jr.loadtimeseries(dataset,'Pressure',3,struc20)
P0 = outdict['clean']
loht_T, hiht_T = jr.interpolationHeights(dataset,ht,'Temperature')
outdict  = jr.loadtimeseries(dataset,'Temperature',loht_T,struc20)
T_lo = outdict['clean']    # celsius
outdict  = jr.loadtimeseries(dataset,'Temperature',hiht_T,struc20)
T_hi = outdict['clean']    # celsius

# interpolate dewpoint temperature in Celsius
DPz_bar = jr.interpolateparameter(dataset,ht,np.nanmean(DP_lo),
                                 np.nanmean(DP_hi),'Dewpt_Temp')
Tz_bar = jr.interpolateparameter(dataset,ht,np.nanmean(T_lo),
                                 np.nanmean(T_hi),'Temperature')
P0_bar = np.nanmean(P0)
DP0_bar = np.nanmean(DP0)
T0_bar_K = np.nanmean(T0) + 273.15
               
# vapor pressure in hPa at height
if (DPz_bar > 0): A, B = 7.5, 237.3
else:            A, B = 9.5, 265.5
e = 6.11 * 10 ** ((DPz_bar*A)/(DPz_bar + B))
   
# vapor pressure in hPa at 3m
if (DP0_bar > 0): A, B = 7.5, 237.3
else:            A, B = 9.5, 265.5
e0 = 6.11 * 10 ** ((DP0_bar*A)/(DP0_bar + B))
               
# extrapolate pressure
g, R = 9.81, 287
q0 = 0.622*e0/P0_bar
Tv0 = T0_bar_K*(1 + 0.61*q0)
dPdz = - (g * P0_bar) / (R * Tv0)
Pz_bar = P0_bar + (ht - 3)*dPdz

# specific humidty at height
qz = 0.622*e/Pz_bar
Tv_K = (Tz_bar + 273.15)*(1 + 0.61*qz)
Tv = Tv_K - 273.15