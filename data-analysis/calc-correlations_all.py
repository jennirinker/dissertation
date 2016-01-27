""" Plot correlation plots
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
import scipy.stats
import scipy.io as scio
import os

# datasets to calculate correlations for
datasets = ['NREL','fluela','PM06','texastech']

# define directory where wind parameters are stored (unused for matlab)
BaseDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                  'processed_data'

# -----------------------------------------------------------------------------

for dataset in datasets:

    # load file with all metadata
    AllHtsName = '{:s}-metadata_allheights.mat'.format(dataset)
    AllHtsPath = os.path.join(BaseDir,AllHtsName)
    outdict = scio.loadmat(AllHtsPath,squeeze_me=True)
    values  = outdict['values']
    IDs = outdict['IDs']
    all_fields  = [s.rstrip() for s in outdict['all_fields']]
    
    # standardize the random variables
    F = (np.arange(values.shape[0])+1.)/(values.shape[0]+1.)
    gvals = np.empty(values.shape)[:,1:]
    for i_parm in range(gvals.shape[1]):
        gvals[np.argsort(values[:,i_parm+1]),i_parm] = scipy.stats.norm.ppf(F)
    
    # calculate Pearson correlation coefficient
    R = np.corrcoef(gvals.T)
    
    # save result
    SaveDict = {}
    SaveDict['IDs']        = IDs
    SaveDict['all_fields'] = all_fields
    SaveDict['R']          = R
    DictName = '{:s}_correlations.mat'.format(dataset)
    DictPath = os.path.join(BaseDir,DictName)
    scio.savemat(DictPath,SaveDict)
    print('Dict {:s} saved'.format(DictName))