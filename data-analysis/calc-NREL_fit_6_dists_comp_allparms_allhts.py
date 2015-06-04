"""
Define list of 6 distributions for wind parameters, fit comp distributions of
those parameters to each parameter at each height, and save in a .txt file.
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
import JR_Library.main as jr

import scipy.io as scio
import numpy as np


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
parms = ['U','sigma_u','L','rho']
Q_Ts = [0.80,0.85,0.90,0.95,1.00]

# initialize wind parameter dist params
p_parms = []

# loop through parameters
for iP in range(4):
    print('Processing parameter {}'.format(parms[iP]))
    
    # initialize height parameters
    h_parms = []
    
    # choose list of distribution candidates
    if (iP < 3): dist_cands = half_cands
    elif (iP == 3): dist_cands = fine_cands    
    
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
            
            # isolate distribution
            dist_name = dist_cands[iD]
            print('      ...{}'.format(dist_name))
                
            # initialize parameter list
            Q_parms = []
        
            # loop through threshold values
            for iQ in range(len(Q_Ts)):
        
                # isolate quantile
                Q_T  = Q_Ts[iQ]
                print('        ...quantile {}'.format(Q_T))
        
                # calculate threshold value
                x, N = np.sort(x), x.size
                if (Q_T == 1.0): x_T = float('inf')
                else:            x_T  = x[N*Q_T]
        
                # optimize distribution parameters
                p_main, p_GP = jr.fitcompositeparameters( \
                    x,dist_name,x_T)
        
                # calculate corresponding NSAE
                NSAE = jr.compositeNSAE(x,dist_name,p_main, \
                                     x_T,p_GP)
        
                # save parameters and NSAE
                fit_parms    = (dist_name,p_main,x_T,p_GP,NSAE)
        
                Q_parms.append(fit_parms)
            
            # make dist list
            d_parms.append(Q_parms)
        
        # save optimum parameters
        h_parms.append(d_parms)
        
    # make full parameter list
    p_parms.append(h_parms)
    
# %% ======================= save data ===============================

choice = input('Do you want to save the distribution information? [1/0]: ')

if choice:
    import json
    fname = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
        'processed_data\\NREL_6dist_comp_parms.txt'
    with open(fname,'w') as f:
        json.dump(p_parms, f)
    print('Parameters saved at:\n{}'.format(fname))
    
# %% ================== calculate/save lowest NSAE ============================

p_parms_opt = []

for iP in range(4):
    
    if (iP < 3): dist_cands = half_cands
    elif (iP == 3): dist_cands = fine_cands 
    
    h_parms_opt = []
    for iH in range(6):
        
        NSAEs = np.empty((len(dist_cands),len(Q_Ts)))
        
        for iD in range(len(dist_cands)):
            for iQ in range(len(Q_Ts)):
                NSAEs[iD,iQ] = p_parms[iP][iH][iD][iQ][4]
                
        iD_min, iQ_min = np.where(NSAEs == NSAEs.min())
        
        opt_parms = p_parms[iP][iH][iD_min][iQ_min]
        
        h_parms_opt.append(opt_parms)
        
    p_parms_opt.append(h_parms_opt)
    
choice = input('Do you want to save the optimal distribution information? [1/0]: ')

if choice:
    import json
    fname = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
        'processed_data\\NREL_optdist_comp_parms.txt'
    with open(fname,'w') as f:
        json.dump(p_parms_opt, f)
    print('Parameters saved at:\n{}'.format(fname))





