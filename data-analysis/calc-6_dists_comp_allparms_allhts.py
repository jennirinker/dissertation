"""
Define list of 6 distributions for wind parameters, fit comp distributions of
those parameters to each parameter at each height, and save in a .txt file.
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
import JR_Library.main as jr

import numpy as np
import os, json


# =========================== user variables ==================================

# choose which dataset
#dataset = 'NREL-mat'
dataset = 'NREL'
#dataset = 'fluela'

# define directory where wind parameters are stored (unused for matlab)
basedir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                  'processed_data'

# parameters to fit distributions to
parms = ['Mean_Wind_Speed','Sigma_u','Tau_u','Concentration_u']

# cutoff threshold values to evaluate
Q_Ts = [0.80,0.85,0.90,0.95,1.00]

# flag to save distribution information and optimal distribution
save_dist = 1

# ============================ load data ======================================

print('\nFitting composite marginal distributions to dataset \"{:s}\"\n'.format(dataset))

if ('mat' in dataset):
    fields, raw_parms = jr.loadNRELmatlab()
    dataset_flag = 'NREL'
else:
    fname = dataset + '-metadata.mat'
    fpath = os.path.join(basedir,fname)
    fields, raw_parms = jr.loadmetadata(fpath)
    dataset_flag = dataset

# screen metadata, get measurement heights and columns for data
clean = jr.screenmetadata(fields,raw_parms,dataset_flag)
heights = jr.datasetSpecs(dataset_flag)['IDs']
htCol  = fields.index('ID')

# ========================= fit distributions =================================

# probability distribution candidates
half_cands = ['lognorm','genextreme',\
    'chi2','gengamma','exponweib','expon']      # U, sigma_u, L
fine_cands = ['anglit','beta', \
    'genextreme','gengamma','lognorm']       # rho

# loop through parameters
p_parms = []
for iP in range(len(parms)):
    parm = parms[iP]
    print('Processing parameter \"{:s}\" ({:d}/{:d})'.format(parm,iP,len(parms)))
    
    # initialize height parameters
    h_parms = []
    
    # choose list of distribution candidates
    if (parm == 'Concentration_u'): dist_cands = fine_cands 
    else:                           dist_cands = half_cands   
    
    # loop through heights
    for iH in range(heights.size):
        ht = heights[iH]
        print('  height {:.1f} m ({:d}/{:d})'.format(ht,iH,heights.size))
                
        # isolate parameters for that height
        idx_ht = np.where(clean[:,htCol]==ht)[0]
        parms_ht = clean[idx_ht,fields.index(parm)]
                
        # sort into increasing order
        x = np.sort(parms_ht)

        # loop through distribution candidates
        d_parms = []
        for iD in range(len(dist_cands)):
            
            # isolate distribution
            dist_name = dist_cands[iD]
            print('      {}'.format(dist_name))
                
            # loop through threshold values
            Q_parms = []
            for iQ in range(len(Q_Ts)):
        
                # isolate quantile
                Q_T  = Q_Ts[iQ]
                print('        quantile {}'.format(Q_T))
        
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

    
# ================== calculate/save lowest NSAE ============================

# loop through parameters
p_parms_opt = []
for iP in range(len(parms)):
    
    # choose list of distribution candidates
    if (parm == 'Concentration_u'): dist_cands = fine_cands 
    else:                           dist_cands = half_cands   
    
    # loop through heights
    h_parms_opt = []
    for iH in range(heights.size):
        
        # create NSAE matrix
        NSAEs = np.empty((len(dist_cands),len(Q_Ts)))
        for iD in range(len(dist_cands)):
            for iQ in range(len(Q_Ts)):
                NSAEs[iD,iQ] = p_parms[iP][iH][iD][iQ][4]
                
        # find distribution/threshold index for smallest NSAE
        iD_min, iQ_min = np.where(NSAEs == NSAEs.min())
        
        # get optimal parameters 
        opt_parms = p_parms[iP][iH][iD_min][iQ_min]
        
        h_parms_opt.append(opt_parms)
        
    p_parms_opt.append(h_parms_opt)
    

print('\nProcessing complete.\n')    
    
# ========================= save distribution info ============================
    
if save_dist:
    
    # create dictionary for probability parameters
    dist_dict                = {}
    dist_dict['parms']       = parms
    dist_dict['qts']         = Q_Ts
    dist_dict['p_parms']     = p_parms
    dist_dict['p_parms_opt'] = p_parms_opt
    
    fname = '{:s}_6dist_comp_parms.txt'.format(dataset)
    fpath = os.path.join(basedir,fname)
    with open(fpath,'w') as f:
        json.dump(dist_dict, f)
    print('Parameters saved at:\n{:s}\n'.format(fpath))





