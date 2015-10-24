"""
Plots of M4 data to figure out which values to simulate in FAST on 
Peregrine
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
import scipy.io as scio
import matplotlib.pyplot as plt
import scipy.stats

# save figure?
saveFig = 0

# load data
fname = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
        'processed_data\\NREL-metadata_allheights.mat'
struc = scio.loadmat(fname)
fields = [s.strip() for s in struc['fields']]
data   = struc['values'][:,1:]
hts    = struc['heights']

# get parameter indices
i_U   = fields.index('U')
i_sig = fields.index('sig_u')
i_L   = fields.index('L_u')
i_rho = fields.index('rho_u')

# get data
U     = data[:,i_U]
sig_u = data[:,i_sig]
I     = sig_u/U
L     = data[:,i_L]
rho   = data[:,i_rho]

# initialize useful things
plt.figure(1,figsize=(6.5,6.5))
plt.clf()
F = (1+np.arange(rho.size))/(1.+rho.size)

# turbulence intensity
plt.subplot(3,2,1)
plt.hist(I,histtype='step',normed=True,bins=20)
plt.title('TI Histogram',fontsize='small')

plt.subplot(3,2,2)
plt.plot(F,np.sort(I))
plt.grid('on')
plt.title('TI ICDF',fontsize='small')

# Kaimal length scale
plt.subplot(3,2,3)
plt.hist(np.log10(L),histtype='step',normed=True,bins=20)
plt.title('Kaimal Length Scale Histogram (log10)',fontsize='small')

plt.subplot(3,2,4)
plt.plot(F,np.sort(np.log10(L)))
plt.grid('on')
plt.title('Kaimal Length Scale ICDF (log10)',fontsize='small')

# get parameter indices
plt.subplot(3,2,5)
plt.hist(rho,histtype='step',normed=True,bins=20)
plt.title('Concentration Parameter Histogram',fontsize='small')

plt.subplot(3,2,6)
plt.plot(F,np.sort(rho))
plt.grid('on')
plt.title('Concentration Parameter ICDF',fontsize='small')

plt.tight_layout()

# save figure if requested
if saveFig:
    plt.savefig('plot-wind_parameter_sampling.png')
    print 'Figure saved.'
    