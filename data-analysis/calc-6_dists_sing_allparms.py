"""
Define list of 6 distributions for wind parameters, fit single distributions of
those parameters to each parameter at each height, and save in a .txt file.

Does NOT require data to be present at all heights
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
#dataset = 'NREL'
dataset = 'fluela'
#dataset = 'PM06'

# define directory where wind parameters are stored (unused for matlab)
BaseDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                  'processed_data'

# parameters to fit distributions to
parms = ['Mean_Wind_Speed','Sigma_u','Tau_u','Concentration_u']

# flag to save distribution information
save_dist = 1

# ============================ load data ======================================

print('\nFitting single marginal distributions to dataset \"{:s}\"\n'.format(dataset))

# load metadata
if ('mat' in dataset):
    fields, raw_parms = jr.loadNRELmatlab()
    dataset_flag = 'NREL'
else:
    fname = dataset + '-metadata.mat'
    fpath = os.path.join(BaseDir,fname)
    fields, raw_parms = jr.loadmetadata(fpath)
    dataset_flag = dataset

# screen metadata, get measurement heights and columns for data
clean = jr.screenmetadata(fields,raw_parms,dataset_flag)
IDs   = jr.datasetSpecs(dataset_flag)['IDs']
htCol  = fields.index('ID')

# ========================= fit distributions ===============================

# probability distribution candidates
half_cands = ['lognorm','genextreme',\
    'chi2','gengamma','exponweib','expon']      # U, sigma_u, L
fine_cands = ['anglit','beta', \
    'genextreme','gengamma','lognorm']       # rho
x_T = float('inf')                              # fitting single distributions

# initialize wind parameter dist params
p_parms = []

# loop through parameters
for iP in range(len(parms)):
    parm = parms[iP]
    print('Processing parameter \"{:s}\" ({:d}/{:d})'.format(parm,iP,len(parms)))
        
    # initialize height parameters
    h_parms = []
    
    # choose list of distribution candidates
    if (parm == 'Concentration_u'): dist_cands = fine_cands 
    else:                           dist_cands = half_cands    
    
    # loop through IDs
    for iH in range(len(IDs)):
        ht = IDs[iH]
        print('  height {:.1f} m ({:d}/{:d})'.format(ht,iH,len(IDs)))
        
        # isolate parameters for that height
        idx_ht = np.where(clean[:,htCol]==ht)[0]
        parms_ht = clean[idx_ht,fields.index(parm)]
                
        # sort into increasing order
        x = np.sort(parms_ht)

        # initialize list of distribution parameters
        d_parms = []
        
        # loop through distribution candidates
        for iD in range(len(dist_cands)):
                
            # separate distribution name
            dist_name = dist_cands[iD]
            print('    {}'.format(dist_name))
            #    
            # fit distribution to data
            p_main, p_GP = jr.fitcompositeparameters(x,dist_name,x_T)
        
            # calculate NSAE
            NSAE = jr.compositeNSAE(x,dist_name,p_main,x_T,p_GP)
            
            # create parameter list
            fit_parms = (dist_name,p_main,x_T,p_GP,NSAE)
            
            # make dist list
            d_parms.append(fit_parms)
        
        # make height list
        h_parms.append(d_parms)
        
    # make full parameter list
    p_parms.append(h_parms)
    
print('\nProcessing complete.\n')
    
# ================== calculate/save lowest NSAE ============================

# loop through parameters
p_parms_opt = []
for iP in range(len(parms)):
    
    # choose list of distribution candidates
    if (parm == 'Concentration_u'): dist_cands = fine_cands 
    else:                           dist_cands = half_cands   
    
    # loop through heights
    h_parms_opt = []
    for iH in range(len(IDs)):
        
        # create NSAE matrix
        NSAEs = np.empty((len(dist_cands)))
        for iD in range(len(dist_cands)):
            NSAEs[iD] = p_parms[iP][iH][iD][4]
                
        # find distribution/threshold index for smallest NSAE
        iD_min = NSAEs.argmin()
        
        # get optimal parameters 
        opt_parms = p_parms[iP][iH][iD_min]
        
        h_parms_opt.append(opt_parms)
        
    p_parms_opt.append(h_parms_opt)
    

print('\nProcessing complete.\n')     
    
#  ======================= save data ===============================

if save_dist:
    
    # create dictionary for probability parameters
    dist_dict                = {}
    dist_dict['parms']       = parms
    dist_dict['p_parms']     = p_parms
    dist_dict['p_parms_opt'] = p_parms_opt
    
    fname = '{:s}_6dist_sing_parms.txt'.format(dataset)
    fpath = os.path.join(BaseDir,fname)
    with open(fpath,'w') as f:
        json.dump(dist_dict, f)
    print('Parameters saved at:\n{:s}\n'.format(fpath))




