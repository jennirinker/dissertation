"""
Plot the errors for different distributions, heights, and parameters in a 
scatter plot to clearly show which distributions generally fit parameters the
best.
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
import JR_Library.main as jr
#
import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt
import json

# %%======================== load data/parameters =============================

# path to matlab-processed metadata table
matpath = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
    '2015-02-28_temporal coherence in data\\code\\tempStruc.mat'

# load the structure
struc = scio.loadmat(matpath)

# extract the information
pdfTable = struc['tempStruc'][0,0][3]

# load fit parameters
fname = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
        'processed_data\\NREL_6dist_parms.txt'
with open(fname,'r') as f:
    p_parms = json.load(f)
    
# %% ======================= plot distributions ===============================

# list of strings of variables
parm_str = ['U','sigma_u','L','rho']

# initialize figure
plt.figure(1,figsize=(15,12))
plt.clf()

# create empirical distribution
N = pdfTable.shape[0]
F_emp = np.arange(N)/(N+1.)

# plot each parameter
iH = 0
for iP in range(4):
    
    # extract data
    x = pdfTable[:,5*iH + iP]
    x = np.sort(x)
    
    # empirical distribution
    plt.subplot(2,2,iP+1)
    plt.plot(x,F_emp)

    # create legend string
    legstr = ['empirical (0.0)']
    for iD in range(len(p_parms[iP][iH])):
        F = jr.compositeCDF(x,p_parms[iP][iH][iD][0],p_parms[iP][iH][iD][1], \
            p_parms[iP][iH][iD][2],p_parms[iP][iH][iD][3])
        plt.plot(x,F)
    
        leg_entry = p_parms[iP][iH][iD][0] + ' ({:.4f})' \
            .format(p_parms[iP][iH][iD][4])
        legstr = legstr + [leg_entry]
    
    plt.legend(legstr,loc='lower right')
    plt.title('Parameter {}'.format(parm_str[iP]))
