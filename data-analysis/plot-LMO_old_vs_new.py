"""
Load entire metadata set, screen it, isolate records with clean data at all
heights, and save the result.
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import os
import numpy as np
import matplotlib.pyplot as plt

# choose which dataset to process
dataset = 'NREL'

# define directory where wind parameters are stored (unused for matlab)
basedir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                  'processed_data'

# import raw metadata
fname = dataset + '-metadata.mat'
fpath = os.path.join(basedir,fname)
fields, raw_parms = jr.loadmetadata(fpath)

# get columns
htCol  = fields.index('Height')
UCol   = fields.index('Mean_Wind_Speed')
sigCol = fields.index('Sigma_u')
tauCol = fields.index('tau_u')
rhoCol = fields.index('Concentration_u')
muCol  = fields.index('Location_u')
LMOoldCol = fields.index('MO_Length_interp')
LMOnewCol = fields.index('MO_Length_virt')

# screen metadata, get measurement heights
clean   = jr.screenmetadata(fields,raw_parms,dataset)
heights = jr.datasetSpecs(dataset)[2]

# get values
LMOold = clean[:,LMOoldCol]
LMOnew = clean[:,LMOnewCol]
zetaold = clean[:,htCol]/LMOold
zetanew = clean[:,htCol]/LMOnew

# plot comparison
Lmax = 5
plt.figure(1,figsize=(4,4))
plt.clf()
plt.plot([-Lmax,Lmax],[-Lmax,Lmax],'k:')
plt.scatter(zetaold,zetanew,s=1)
plt.xlim([-Lmax,Lmax])
plt.ylim([-Lmax,Lmax])