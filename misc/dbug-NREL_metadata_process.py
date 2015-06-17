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

# %% ============================= load data ==================================

## load fields, metadata from matlab
#flds_mat, md_mat = jr.loadNRELmatlab()
#
## load python-processed metadata
#fname = 'C:\\Users\\jrinker\\Dropbox\\research' + \
#    '\\processed_data\\NREL_metadata.mat'
#struc = scio.loadmat(fname)
#flds_py = struc['fields']
#md_py   = struc['metadata']

# %% ========================== compare metadata ==============================

# convert first python date to datevec
idx_py = 3
time_flt = md_py[idx_py,0]
time_tup = jr.timeflt2tup(time_flt)
ht       = int(md_py[idx_py,2])
rec_vec  = np.asarray(time_tup + (ht,))
print(rec_vec)

# find index of that date in matlab metadata
idx_mat = int(np.squeeze(np.where(np.all(md_mat[:,:6]==rec_vec,axis=1))))

# get vectors of data
dat_py  = md_py[0]
dat_mat = md_mat[idx_mat,:]

# manually extract parameters
row = jr.extractNRELparameters(struc,ht)
row[0] = time_flt
row[1] = calendar.timegm(time.gmtime())
row[2] = ht

# compare parameters
parms = ['WS_Cup','Dir','Prec','U','sig_u','rho_u','mu_u','sig_v','rho_v', \
    'mu_v','sig_w','rho_w','mu_w','tau_u','tau_v','tau_w']
parms_py  = np.append(dat_py[3:16],dat_py[20:])
parms_py2 = np.append(row[3:16],row[20:])
parms_mat = dat_mat[[6,7,8,13,14,15,16,17,18,19,20,21,22,27,28,29]]



#
print('    '.join(parms))
print('----------------------------------------------------------------------')
print('   '.join(['{:.3f}'.format(i) for i in parms_py]))
print('   '.join(['{:.3f}'.format(i) for i in parms_mat]))
print('   '.join(['{:.3f}'.format(i) for i in parms_py2]))
print('')
print('   '.join(['{:.3f}'.format(i) for i in parms_mat-parms_py2]))

## manually extract parameters
#t, u, v, w = jr.loadtimeseries('NREL',time_flt,ht)
#uhat = scipy.signal.detrend(u) + np.mean(u)
#vhat = scipy.signal.detrend(v)
#what = scipy.signal.detrend(w)
#
#plt.figure(1)
#plt.clf()
#plt.subplot(2,1,1)
#plt.plot(t,u)
#plt.plot(t,v)
#plt.plot(t,w)
#plt.plot(t,uhat)
#plt.legend(['u','v','w','uhat'])
#plt.subplot(2,1,2)
#plt.plot(t,uhat*vhat)
#plt.plot(t,uhat*what)
#plt.plot(t,vhat*what)
#plt.legend(['uv','uw','vw'])

## set first date from matlab, find index in matlab array
#day1 = [2012,02,13,16,30]
#day_idx = np.where(np.all(md_mat[:,:5]==[2012,02,13,16,30],axis=1))[0][0]
#row_mat = md_mat[day_idx,13:]
#
## load time history
#fpath = 'G:\\data\\nrel-20Hz\\2012\\02\\13\\02_13_2012_16_30_00_038.mat'
#struc = scio.loadmat(fpath)
#row_py = jr.extractNRELparameters(struc,15)[:,6:]
#flds_py = jr.metadataFields('NREL')
#
## print percent differences in u:
#print('Percent differences in u, sigma_u, rho, and mu:')
#print(np.squeeze((row_mat[:4]-row_py[0,:4])/row_mat[:4]*100))

