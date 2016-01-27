"""
Summarize processed wind parameters from a given dataset
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as scio

# choose plot style
plt.style.use('C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
            '2016-02-15_dissertation\\figure_code\\duke_paper.mplstyle')

# choose which dataset
dataset, fignum = 'texastech', 1

# define directory where wind parameters are stored (unused for matlab)
basedir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                  'processed_data'
                  
print('\nProcessing dataset {}'.format(dataset))

if ('mat' in dataset):
    fields, raw_parms = jr.loadNRELmatlab()
    dataset = 'NREL'
else:
#    fname = dataset + '-metadata.mat'
#    fpath = os.path.join(basedir,fname)
    fields, clean = jr.loadmetadata(dataset)

print(np.nanmin(clean[:,fields.index('ID')]))

# screen metadata
screen = jr.screenmetadata(fields,clean,dataset)

print(screen[:,fields.index('ID')].min())
#    
## get columns
#UsCol   = fields.index('Wind_Speed_Sonic')
#UpCol   = fields.index('Wind_Speed_UVW')
#siguCol = fields.index('Sigma_u')
#rhouCol = fields.index('Concentration_u')
    
# remove NaN values from wind speeds
#clean = jr.cleantexastech(raw_parms,fields)

## define parameters for simplicity
#Us   = raw_parms[:,UsCol]
#Up   = raw_parms[:,UpCol]
#sigu = raw_parms[:,siguCol]
#rhou = raw_parms[:,rhouCol]
#
#plt.figure(1)
#plt.clf()
#
#plt.subplot(2,2,1)
#plt.hist(Us,bins=20)
#plt.title('Sonic WS')
#
#plt.subplot(2,2,2)
#plt.hist(Up,bins=20)
#plt.title('UVW WS')
#
#plt.subplot(2,2,3)
#plt.hist(sigu,bins=20)
#plt.title('Sigma_u')
#
#plt.subplot(2,2,4)
#plt.hist(rhou,bins=20)
#plt.title('Rho_u')
#
#plt.tight_layout()

## plot data with highest rho
#idx = np.where(clean[:,UsCol] == np.sort(Us)[-1])[0][0]
#ID,time_flt = clean[idx,fields.index('ID')],clean[idx,0]
#x = jr.loadtimeseries(dataset,'Sonic_x',ID,time_flt)['clean']
#y = jr.loadtimeseries(dataset,'Sonic_y',ID,time_flt)['clean']
#z = jr.loadtimeseries(dataset,'Sonic_z',ID,time_flt)['clean']
#print(jr.loadtimeseries(dataset,'Sonic_x',ID,time_flt)['flags'])
#plt.figure(2)
#plt.clf()
#plt.plot(x)
#plt.plot(y)
#plt.plot(z)
