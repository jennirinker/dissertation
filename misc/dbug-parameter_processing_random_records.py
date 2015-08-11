"""
Testing the parallel wind parameter processing routine in Python
"""
import sys, os
from joblib import Parallel, delayed
import json
import numpy as np
import scipy.io as scio

# clean date 1: (2013,3,11,6,20)
# clean date 2: (2013,4,15,3,0)
# clean date 3: (2013,4,15,3,10)
# clean date 4: (2013,5,30,3,50)
# clean date 5: (2013,11,3,14,30)

# =================== COMPARE MAT/PY VALUES FOR CLEAN ======================

## calculate python metadata
#i_files = []
#timestamp = (2013,3,11,6,20)
#fpath     = jr.NRELtime2fpath(timestamp)
#i_files.append(list_mats.index(fpath))
#md_list = []
#for i in range(len(i_files)):
#    md_list.append(jr.NRELlistmetadata(i_files[i],list_mats))
#md_py = np.empty((6*len(i_files),len(fields)))
#for i in range(len(i_files)):
#    md_py[6*i:6*(i+1)] = np.asarray(md_list[i])
#dat_py = md_py[0,:]
#
## get matlab index/metadata
#flds_mat, raw_mat = jr.loadNRELmatlab()
#clean_mat = jr.screenmetadata(flds_mat,raw_mat,'NREL')
#ht = 15
#rec_vec  = np.asarray(timestamp + (ht,))         # time/height tuple
#idx_mat = np.squeeze(np.where(np.all( \
#    clean_mat[:,:6]==rec_vec,axis=1)))             # corresponding matlab index
#dat_mat = clean_mat[idx_mat,:]
#
#
#
#parms = ['WS_Cup','Dir  ','Prec ','U   ','sig_u ','rho_u','mu_u','sig_v','rho_v', \
#    'mu_v','sig_w','rho_w','mu_w','tau_u','tau_v','tau_w']
#prms_mdpy  = np.append(dat_py[3:16],dat_py[20:])
#prms_mdmat = dat_mat[[6,7,8,13,14,15,16,17,18,19,20,21,22,27,28,29]]
#
#
#print('   '.join(parms))
#print('----------------------------------------------------------------------')
#print('   '.join(['{:.3f}'.format(i) for i in prms_mdmat]))
#print('   '.join(['{:.3f}'.format(i) for i in prms_mdpy]))


# ================= CALCULATE A FEW METADATAS IN PARALLEL ====================

if (__name__ == '__main__'):
    libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
    if (libpath not in sys.path): sys.path.append(libpath)
    import JR_Library.main as jr
    
    # load list of mat files
    basedir  = jr.getBasedir('NREL')
    fields   = jr.metadataFields('NREL')
    lmats_fname = [fp for fp in os.listdir(basedir) if 'listmats' in fp][0]
    lmats_fpath = os.path.join(basedir,lmats_fname)
    with open(lmats_fpath,'r') as f:
        list_mats = json.load(f)

    
    # get list of file indices
    i_files = []
    timestamp = (2013,3,11,6,20)
    fpath     = jr.NRELtime2fpath(timestamp)
    i_files.append(list_mats.index(fpath))
    timestamp = (2013,11,3,14,30)
    fpath     = jr.NRELtime2fpath(timestamp)
    i_files.append(list_mats.index(fpath))
    
    # variables for later
    basedir  = jr.getBasedir('NREL')
    njobs    = min(4,len(i_files))
    fields   = jr.metadataFields('NREL')
        
    # process files in parallel
    md_list = []
#    for i in range(len(i_files)):
#        md_list.append(jr.NRELlistmetadata(i_files[i],list_mats))
    md_list = Parallel(n_jobs=njobs,verbose=9) \
        (delayed(jr.NRELlistmetadata)(i_files[i],list_mats) for i \
            in range(len(i_files)))
        
    # convert list of arrays to metadata
    metadata = np.empty((6*len(i_files),len(fields)))
    for i in range(len(i_files)):
        metadata[6*i:6*(i+1)] = np.asarray(md_list[i])




