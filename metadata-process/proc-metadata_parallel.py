"""
Script to process the NREL metadata in parallel.
"""

import sys, os
from joblib import Parallel, delayed
import json
import numpy as np
import scipy.io as scio

if (__name__ == '__main__'):
    libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
#    libpath = 'E:\\NREL-metadata-process'
    if (libpath not in sys.path): sys.path.append(libpath)
    import JR_Library.main as jr
    
    # variables for later
#    datasets       = ['NREL','fluela','PM06','texastech']
    dataset        = 'texastech'
    basedir        = jr.getBasedir(dataset)
    njobs          = 2
    fields         = jr.metadataFields(dataset)
    save_data      = 0
    in_parallel    = 0
    
    # see if list of mat files exists, create it if it doesn't
    lmats_fname = [fp for fp in os.listdir(basedir) if 'listmats' in fp]
    if not lmats_fname:
        jr.list_matfiles(basedir,save=1)
        lmats_fname = [fp for fp in os.listdir(basedir) if 'listmats' in fp]
    lmats_fname = lmats_fname[0]
        
    # load list of mat files
    lmats_fpath = os.path.join(basedir,lmats_fname)
    with open(lmats_fpath,'r') as f:
        list_mats = json.load(f)
    n_files = len(list_mats)
    
    # number of files to process (make small for debugging)
    n_proc = n_files                            # no. files to process
    n_proc = 2                            # no. files to process
        
    # process files in parallel
    if in_parallel:
        md_list = Parallel(n_jobs=njobs,verbose=9) \
            (delayed(jr.listmetadata)(dataset,i,list_mats) for i in range(n_proc))

    # process files NOT in parallel (mostly for debugging)
    else:        
        md_list = []
        for i in range(n_proc):
            tmp = jr.listmetadata(dataset,i,list_mats)
            md_list.append(tmp)
        
    # convert list of arrays to metadata
    n_IDs = len(jr.datasetSpecs(dataset)['IDs'])
    metadata = np.empty((n_IDs*n_proc,len(fields)))
    for i in range(n_proc):
        metadata[n_IDs*i:n_IDs*(i+1)] = np.asarray(md_list[i])
        
    # save output
    if save_data:
        foutname = dataset + '-metadata.mat'
        foutpath = os.path.join('C:\\Users\\jrinker\\Dropbox\\' + \
                                    'research\\processed_data',foutname)
        mdict             = {}
        mdict['fields']   = fields
        mdict['metadata'] = metadata
        scio.savemat(foutpath,mdict)
        
        print('\nData saved to {:s}'.format(foutpath))
        
    print('\nScript complete.\n')
