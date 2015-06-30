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
#flds_mat, md_mat = jr.loadNRELmatlab()
#
## load python-processed metadata
#fname = 'C:\\Users\\jrinker\\Dropbox\\research' + \
#    '\\processed_data\\NREL_metadata.mat'
#struc = scio.loadmat(fname)
#flds_py = struc['fields']
#md_py   = struc['metadata']

# %% ====================== compare metadata values ===========================

# get python and matlab indices
#idx_rand = 1                                      # python index
idx_rand = metadata[:,8].argmax()
time_flt = metadata[idx_rand,0]                   # float of timestamp
time_tup = jr.timeflt2tup(time_flt)             # tuple of timestamp
ht       = int(metadata[idx_rand,2])             # measurment height
rec_vec  = np.asarray(time_tup + (ht,))         # time/height tuple
idx_mat = np.squeeze(np.where(np.all( \
    md_mat[:,:6]==rec_vec,axis=1)))             # corresponding matlab index
idx_py = np.squeeze(np.where(np.logical_and( \
    time_flt == md_py[:,0],ht == md_py[:,2])))     # corresponding matlab index

# if the record exists in the matlab metadata
if (idx_mat.size > 0):
    
    # get values from the processed metadata arrays
    dat_py  = md_py[idx_py,:]
    dat_mat = md_mat[idx_mat,:]
    dat_py2 = metadata[idx_rand,:]
    
    # print time stamps to ensure correct indices
    print('rand_py date',jr.timeflt2tup(dat_py2[0]))
    print('md_py date',jr.timeflt2tup(dat_py[0]))
    print('md_mat date',dat_mat[:6])

    # rearrange parameters from different metadata
    parms = ['WS_Cup','Dir  ','Prec ','U   ','sig_u ','rho_u','mu_u','sig_v','rho_v', \
        'mu_v','sig_w','rho_w','mu_w','tau_u','tau_v','tau_w']
    prms_mdpy  = np.append(dat_py[3:16],dat_py[20:])
    prms_mdmat = dat_mat[[6,7,8,13,14,15,16,17,18,19,20,21,22,27,28,29]]
    
    # rearrange parameters from manual calculations
    prms_mdpy2 = np.append(dat_py2[3:16],dat_py2[20:-2])
    
    # extract parameters manually
    fpath20 = jr.NRELtime2fpath(time_flt)
    struc20 = scio.loadmat(fpath20)
    row     = jr.extractNRELparameters(struc20,ht)
    prms_manpy = np.append(row[3:16],row[20:-2])


    print('   '.join(parms))
    print('----------------------------------------------------------------------')
    print('   '.join(['{:.3f}'.format(i) for i in prms_mdmat]))
    print('   '.join(['{:.3f}'.format(i) for i in prms_mdpy]))
    print('   '.join(['{:.3f}'.format(i) for i in prms_mdpy2]))
    print('   '.join(['{:.3f}'.format(i) for i in prms_manpy]))
#    print('')
#    print('   '.join(['{:.3f}'.format(i) for i in prms_mdpy2-prms_mdmat]))


t, uraw, vraw, wraw = jr.loadtimeseries('NREL',time_flt,ht)
wind_dir = struc20['Cup_WS_10m'][0,0][0]
plt.figure(1)
plt.clf()
plt.plot(t,wind_dir)
plt.title('Cup WS (10 m), 2013/01/10 23:10')
u = jr.cleantimeseries(t,uraw)
plt.figure(2)
plt.clf()
plt.plot(t,u)

## get python and matlab indices
#idx_py = 10                                      # python index
#time_flt = md_py[idx_py,0]                      # float of timestamp
#time_tup = jr.timeflt2tup(time_flt)             # tuple of timestamp
#ht       = int(md_py[idx_py,2])                 # measurment height
#rec_vec  = np.asarray(time_tup + (ht,))         # time/height tuple
#print(rec_vec)
#idx_mat = np.squeeze(np.where(np.all( \
#    md_mat[:,:6]==rec_vec,axis=1)))             # corresponding matlab index
#
## if the record exists in the matlab metadata
#if (idx_mat.size > 0):
#    
#    # get values from the processed metadata arrays
#    dat_py  = md_py[0]
#    dat_mat = md_mat[idx_mat,:]
#
#    # extract parameters manually
#    fpath20 = jr.NRELtime2fpath(time_flt)
#    struc20 = scio.loadmat(fpath20)
#    row     = jr.extractNRELparameters(struc20,ht)
#    row[0]  = time_flt
#    row[1]  = calendar.timegm(time.gmtime())
#    row[2]  = ht
#
#    # rearrange parameters from different metadata
#    parms = ['WS_Cup','Dir  ','Prec ','U   ','sig_u ','rho_u','mu_u','sig_v','rho_v', \
#        'mu_v','sig_w','rho_w','mu_w','tau_u','tau_v','tau_w']
#    prms_mdpy  = np.append(dat_py[3:16],dat_py[20:])
#    prms_mdmat = dat_mat[[6,7,8,13,14,15,16,17,18,19,20,21,22,27,28,29]]
#    
#    # rearrange parameters from manual calculations
#    prms_manpy = np.append(row[3:16],row[20:-2])
#
#
#    print('   '.join(parms))
#    print('----------------------------------------------------------------------')
#    print('   '.join(['{:.3f}'.format(i) for i in prms_mdmat]))
#    print('   '.join(['{:.3f}'.format(i) for i in prms_mdpy]))
#    print('   '.join(['{:.3f}'.format(i) for i in prms_manpy]))
#    print('')
#    print('   '.join(['{:.3f}'.format(i) for i in prms_mdmat-prms_manpy]))


## check wind direction calculations
#WD_26 = struc20['Vane_WD_26m'][0,0][0]
#WD_88  = struc20['Vane_WD_88m'][0,0][0]
#WD_26_mean = np.angle(np.sum(np.exp(1j*WD_26*np.pi/180.)),deg=1)
#WD_88_mean = np.angle(np.sum(np.exp(1j*WD_88*np.pi/180.)),deg=1)
#print(WD_26_mean)
#print(WD_88_mean)


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

