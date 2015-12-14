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

# choose which dataset
#dataset, fignum = 'NREL-mat', 1
#dataset, fignum = 'NREL', 2
#dataset, fignum = 'fluela', 3
dataset, fignum = 'PM06', 4

# define directory where wind parameters are stored (unused for matlab)
basedir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                  'processed_data'
                  
print('\nProcessing dataset {}'.format(dataset))

if ('mat' in dataset):
    fields, raw_parms = jr.loadNRELmatlab()
    dataset = 'NREL'
else:
    fname = dataset + '-metadata.mat'
    fpath = os.path.join(basedir,fname)
    fields, raw_parms = jr.loadmetadata(fpath)

# screen metadata, get measurement heights
clean = jr.screenmetadata(fields,raw_parms,dataset)
heights = jr.datasetSpecs(dataset)['IDs']

# print sizes
n_raw, n_clean = raw_parms.shape[0], clean.shape[0]
print('\nOverview')
print('------------------------------')
print('  Size of raw dataset:   {}'.format(n_raw))
print('  Size of clean dataset: {}'.format(n_clean))
print('  Percent clean:         {:.2f}%'.format(n_clean/float(n_raw)*100.))

# plot distributions
htCol  = fields.index('ID')
#htCol  = fields.index('Height')
UCol   = fields.index('Mean_Wind_Speed')
sigCol = fields.index('Sigma_u')
tauCol = fields.index('Tau_u')
#tauCol = fields.index('tau_u')
rhoCol = fields.index('Concentration_u')
muCol  = fields.index('Location_u')

# initialize figure
plt.figure(fignum,figsize=(6.5,6))
plt.clf()
ax11 = plt.axes([0.12,0.57,0.35,0.35])
ax12 = plt.axes([0.62,0.57,0.35,0.35])
ax21 = plt.axes([0.12,0.12,0.35,0.35])
ax22 = plt.axes([0.62,0.12,0.35,0.35])

print('\nHealthy records by height')
print(  '-------------------------')

for iH in range(len(heights)):
    
    # isolate parameters for that height
    ht = heights[iH]
    idx_ht = np.where(clean[:,htCol]==ht)[0]
    parms_ht = clean[idx_ht,:]
    
    # calculate empirical CDF
    n_recs = parms_ht.shape[0]
    F = np.arange(1,n_recs+1)/(n_recs+1.)
    print('  {:<3}     {}   ({:.1f}%)'.format(ht,n_recs, \
                    (n_recs/float(n_clean)*100)))
    
    # get parameters
    U   = parms_ht[:,UCol]
    sig = parms_ht[:,sigCol]
    tau = parms_ht[:,tauCol]
    L   = tau*U
    rho = parms_ht[:,rhoCol]
    mu  = parms_ht[:,muCol]
    
    # plot CDFs
    ax11.plot(np.sort(U),F,label='{:.0f} m'.format(ht))
    ax12.plot(np.sort(sig),F)
    ax21.plot(np.sort(tau),F)
#    ax21.plot(np.sort(L),F)
    ax22.plot(np.sort(rho),F)
    
# make things pretty
ax21.set_xscale('log')
ax11.set_title('$U$')
ax12.set_title('$\sigma_u$')
ax21.set_title('$L/U$')
#ax21.set_title('$L$')
ax22.set_title(r'$\rho$')
#ax11.legend(loc=4,fontsize='small')
plt.suptitle('Dataset: ' + dataset,fontsize='large')