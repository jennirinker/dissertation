"""
Define list of 6 distributions for wind parameters, fit single distributions of
those parameters to each parameter at each height, and save in a .mat file.
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
import JR_Library.main as jr

import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt


# %%============================= load data ===================================

# path to matlab-processed metadata table
matpath = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
    '2015-02-28_temporal coherence in data\\code\\tempStruc.mat'

# load the structure
struc = scio.loadmat(matpath)

# extract the information
pdfTable = struc['tempStruc'][0,0][3]

# %%========================= fit distributions ===============================

# probability distribution candidates
half_cands = ['lognorm','genextreme',\
    'chi2','gengamma','exponweib','expon']      # U, sigma_u, L
fine_cands = ['anglit','beta', \
    'genextreme','gengamma','lognorm']       # rho
x_T = float('inf')                              # fitting single distributions
parms = ['U','sigma_u','L','rho']

# initialize wind parameter dist params
p_parms = []

# loop through parameters
for iP in range(4):
    print('Processing parameter {}'.format(parms[iP]))
    
    # initialize height parameters
    h_parms = []
    
    # choose list of distribution candidates
    if (iP < 3): dist_cands = half_cands
    elif (iP ==3): dist_cands = fine_cands    
    
    # loop through heights
    for iH in range(6):
        print('  ...height {}'.format(iH))
                
        # extract data
        x = pdfTable[:,5*iH + iP]
        x = np.sort(x)

        # initialize list of distribution parameters
        d_parms = []
        
        # loop through distribution candidates
        for iD in range(len(dist_cands)):
                
            # separate distribution name
            dist_name = dist_cands[iD]
            print('    ...{}'.format(dist_name))
            #    
            # fit distribution to data
            p_main, p_GP = jr.fitcompositeparameters(x,dist_name,x_T)
        
            # calculate NSAE
            NSAE = jr.compositeNSAE(x,dist_name,p_main,x_T,p_GP)
            
            # create parameter list
            fit_parms = [dist_name,p_main,x_T,p_GP,NSAE]
            
            # make dist list
            d_parms.append(fit_parms)
        
        # make height list
        h_parms.append(d_parms)
        
    # make full parameter list
    p_parms.append(h_parms)
    
# %% ======================= save data ===============================

choice = input('Do you want to save the distribution information? [1/0]: ')

if choice:
    import json
    fname = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
        'processed_data\\NREL_6dist_parms.txt'
    with open(fname,'w') as f:
        json.dump(p_parms, f)
    print('Parameters saved at:\n{}'.format(fname))


