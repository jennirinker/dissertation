"""
Debug NREL processing routine in Python
"""
import scipy.io as scio
import sys
sys.path.append('..')
import JR_Library.main as jr
import numpy as np

# load fields, metadata from matlab
flds_mat, md_mat = jr.loadNRELmatlab()

# set first date from matlab, find index in matlab array
day1 = [2012,02,13,16,30]
day_idx = np.where(np.all(md_mat[:,:5]==[2012,02,13,16,30],axis=1))[0][0]
row_mat = md_mat[day_idx,13:]

# load time history
fpath = 'G:\\data\\nrel-20Hz\\2012\\02\\13\\02_13_2012_16_30_00_038.mat'
struc = scio.loadmat(fpath)
row_py = jr.extractNRELparameters(struc,15)[:,6:]
flds_py = jr.metadataFields('NREL')

# print percent differences in u:
print('Percent differences in u, sigma_u, rho, and mu:')
print(np.squeeze((row_mat[:4]-row_py[0,:4])/row_mat[:4]*100))

