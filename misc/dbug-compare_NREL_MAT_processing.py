"""
Debug NREL processing routine in Python
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
import scipy.io as scio
import scipy.signal
import calendar, time

# %% ===================== load metadata structures ===========================

## load fields, metadata from matlab
#flds_mat, raw_mat = jr.loadNRELmatlab()
#clean_mat = jr.screenmetadata(flds_mat,raw_mat,'NREL')
#
## load python-processed metadata
fname = 'C:\\Users\\jrinker\\Dropbox\\research' + \
    '\\processed_data\\NREL-metadata.mat'
flds_py,raw_py = jr.loadmetadata(fname)
clean_py =  jr.screenmetadata(flds_py,raw_py,'NREL')
#del raw_mat, raw_py

# %% ====================== compare metadata values ===========================

# get python and matlab indices
#idx_py = 1                                      # python index
idx_py = clean_py[:,
                  flds_py.index('tau_u')].argmin()
time_flt = clean_py[idx_py,0]                   # float of timestamp
time_tup = jr.timeflt2tup(time_flt)             # tuple of timestamp
ht       = int(clean_py[idx_py,2])             # measurment height
rec_vec  = np.asarray(time_tup + (ht,))         # time/height tuple
idx_mat = np.squeeze(np.where(np.all( \
    clean_mat[:,:6]==rec_vec,axis=1)))             # corresponding matlab index
#idx_py = np.squeeze(np.where(np.logical_and( \
#    time_flt == clean_py[:,0],ht == clean_py[:,2])))     # corresponding matlab index


# if the record exists in the matlab metadata
if (idx_mat.size > 0):
    
    # get values from the processed metadata arrays
    dat_py  = clean_py[idx_py,:]
    dat_mat = clean_mat[idx_mat,:]
    
    # print time stamps to ensure correct indices
    print('clean_py date',time_tup)
    print('clean_mat date',dat_mat[:6])

    # rearrange parameters from different metadata
    parms = ['WS_Cup','Dir  ','Prec ','U   ','sig_u ','rho_u','mu_u','sig_v','rho_v', \
        'mu_v','sig_w','rho_w','mu_w','tau_u','tau_v','tau_w']
    prms_mdpy  = np.append(dat_py[3:16],dat_py[20:])
    prms_mdmat = dat_mat[[6,7,8,13,14,15,16,17,18,19,20,21,22,27,28,29]]
    
    # rearrange parameters from manual calculations
#    prms_mdpy2 = np.append(dat_py2[3:16],dat_py2[20:-2])
    
    # extract parameters manually
    fpath20 = jr.time2fpath('NREL',time_flt)
    struc20 = scio.loadmat(fpath20)
    row     = jr.struc2metadata('NREL',struc20,ht)
    prms_manpy = np.append(row[3:16],row[20:-2])


    print('   '.join(parms))
    print('----------------------------------------------------------------------')
    print('   '.join(['{:.3f}'.format(i) for i in prms_mdmat]))
    print('   '.join(['{:.3f}'.format(i) for i in prms_mdpy]))
    print('   '.join(['{:.3f}'.format(i) for i in prms_manpy]))


else:
    print('Record does not exist in Matlab metadata')