"""
Demo - Andy L_u versus mine
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
import scipy.io as scio
import matplotlib.pyplot as plt

# flag to save figure
saveFig = 0

# load fields, metadata from matlab
flds_mat, raw_mat = jr.loadNRELmatlab()
clean_mat = jr.screenmetadata(flds_mat,raw_mat,'NREL')

# get indices for parameters
Lu_Andy_idx = flds_mat.index('L_u_Andy')
tau_u_idx   = flds_mat.index('tau_u')
U_idx       = flds_mat.index('Mean_Wind_Speed')

# get values necessary for comparison
U        = clean_mat[:,U_idx]
Lu_Andy  = clean_mat[:,Lu_Andy_idx]
Lu_Jenni = clean_mat[:,tau_u_idx]*U

# plot Kaimal length scales
plt.figure(1,figsize=(5,5))
plt.clf()

#plt.hist2d(Lu_Andy,Lu_Jenni,bins=50)

plt.scatter(Lu_Andy,Lu_Jenni,c=U,s=8,edgecolor='none')
plt.xlim([0,1e4])
plt.ylim([0,1e4])
plt.xlabel('Kaimal Length Scale - Andy')
plt.ylabel('Kaimal Length Scale - Jenni')
plt.colorbar()
plt.tight_layout()

# save figure if requested
if saveFig:
    plt.savefig('demo-Andy_vs_Jenni_Lu.png')
    print 'Figure saved.'