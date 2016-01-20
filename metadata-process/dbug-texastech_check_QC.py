"""
debugging texas tech metadata processing
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import scipy.io as scio
import os
import numpy as np
import matplotlib.pyplot as plt

# define dataset
dataset = 'texastech'

BaseDataDir = '.'

# -----------------------------------------------------------------------------

# load metadata
DataName = '{:s}-metadata.mat'.format(dataset)
DataPath = os.path.join(BaseDataDir,DataName)
fields,metadata = jr.loadmetadata(DataPath)

# screen data by sonic speed and wind direction
cleandata = jr.screenmetadata(fields,metadata,dataset)

# extract largest wind speed/sig/etc
U = cleandata[:,fields.index('Mean_Wind_Speed')]
sig = cleandata[:,fields.index('Sigma_u')]
rho = cleandata[:,fields.index('Concentration_u')]

# get time series corresponding to max value(s)
idx = np.where(rho == np.sort(rho)[-1])[0][0]
print(rho[idx])
ID = cleandata[idx,fields.index('ID')]

fpath = jr.time2fpath(dataset,cleandata[idx,fields.index('Record_Time')])
struc_hf = scio.loadmat(fpath,squeeze_me=True)

parms = jr.struc2metadata(dataset,struc_hf,ID)

RecTime = cleandata[idx,fields.index('Record_Time')]
TSfield = 'Sonic_x'
TSdict = jr.loadtimeseries(dataset,TSfield,ID,RecTime)
t, x = TSdict['time'],TSdict['clean']
print(idx,TSdict['flags'])
TSfield = 'Sonic_y'
TSdict = jr.loadtimeseries(dataset,TSfield,ID,RecTime)
t, y = TSdict['time'],TSdict['clean']
print(idx,TSdict['flags'])
TSfield = 'Sonic_z'
TSdict = jr.loadtimeseries(dataset,TSfield,ID,RecTime)
t, z = TSdict['time'],TSdict['clean']
print(idx,TSdict['flags'])

plt.clf()
plt.plot(t,x)
plt.plot(t,y)
plt.plot(t,z)
#
#x_cl, n_spikes = jr.cleantimeseries(t,TSdict['raw'])

