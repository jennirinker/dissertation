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
    dataset  = 'NREL'
    basedir  = jr.getBasedir(dataset)
    njobs    = 4
    fields   = jr.metadataFields(dataset)
    foutname = dataset + '-metadata.mat'
    
    # load saved list of mat files
    lmats_fname = [fp for fp in os.listdir(basedir) if 'listmats' in fp][0]
    lmats_fpath = os.path.join(basedir,lmats_fname)
    with open(lmats_fpath,'r') as f:
        list_mats = json.load(f)
    n_files = len(list_mats)
        
    # process files in parallel
#    md_list = Parallel(n_jobs=njobs,verbose=9) \
#        (delayed(jr.NRELlistmetadata)(i,list_mats) for i in range(n_files))
    n_proc = 1
    md_list = Parallel(n_jobs=njobs,verbose=9) \
        (delayed(jr.listmetadata)(dataset,i,list_mats) for i in range(n_proc))
        
    # convert list of arrays to metadata
    n_heights = jr.datasetSpecs(dataset)[2].size
    metadata = np.empty((n_heights*n_proc,len(fields)))
    for i in range(n_proc):
        metadata[n_heights*i:n_heights*(i+1)] = np.asarray(md_list[i])
        
    # save output
#    mdict             = {}
#    mdict['fields']   = fields
#    mdict['metadata'] = metadata
#    scio.savemat(foutname,mdict)
