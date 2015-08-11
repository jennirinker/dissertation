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
    basedir  = jr.getBasedir('NREL')
    njobs    = 60
    fields   = jr.metadataFields('NREL')
    foutname = 'NREL-metadata.mat'
    
    # load saved list of mat files
    lmats_fname = [fp for fp in os.listdir(basedir) if 'listmats' in fp][0]
    lmats_fpath = os.path.join(basedir,lmats_fname)
    with open(lmats_fpath,'r') as f:
        list_mats = json.load(f)
    n_files = len(list_mats)
        
    # process files in parallel
    md_list = Parallel(n_jobs=njobs,verbose=9) \
        (delayed(jr.NRELlistmetadata)(i,list_mats) for i in range(n_files))
        
    # convert list of arrays to metadata
    metadata = np.empty((6*n_files,len(fields)))
    for i in range(n_files):
        metadata[6*i:6*(i+1)] = np.asarray(md_list[i])
        
    # save output
    mdict             = {}
    mdict['fields']   = fields
    mdict['metadata'] = metadata
    scio.savemat(foutname,mdict)
