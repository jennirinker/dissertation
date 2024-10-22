"""
Compare the empirical CDFs and the best fit CDFs
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
import JR_Library.main as jr

import numpy as np
import os, json
import matplotlib.pyplot as plt

# choose which dataset
datasets = ['NREL','fluela','PM06','texastech']

# define directory where wind parameters are stored (unused for matlab)
basedir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                  'processed_data'

# flag to plot single or composite distribution
#dist_type = 'sing'
dist_type = 'comp'

# plot style
plt.style.use(jr.stylepath('duke_paper'))
FigNum = 1

# -----------------------------------------------------------------------------

iPlot = 0
for dataset in datasets:

    print('\nPlotting {:s} marginal '.format(dist_type) + \
            'distributions to dataset \"{:s}\"\n'.format(dataset))
    
    if ('mat' in dataset):
        fields, clean = jr.loadNRELmatlab()
        dataset_flag = 'NREL'
    else:
        fields, clean = jr.loadmetadata(dataset)
        dataset_flag = dataset
    
    # screen metadata, get measurement heights and columns for data
    screen = jr.screenmetadata(fields,clean,dataset_flag)
    heights = jr.datasetSpecs(dataset_flag)['IDs']
    htCol  = fields.index('ID')
    
    # load distribution information
    dist_fname = '{:s}_6dist_{:s}_parms.txt'.format(dataset,dist_type)
    dist_fpath = os.path.join(basedir,dist_fname)
    with open(dist_fpath,'r') as f:
        dist_dict  = json.load(f)
    p_parms_opt = dist_dict['p_parms_opt']
    parms       = dist_dict['parms']
    
    parms[2] = 'Tau_u'
    
    # ========================== plot distributions ===============================
    
    # initialize figure
    plt.figure(FigNum + iPlot,figsize=(6.,5.5))
    plt.clf()
    ax = []
    ax.append(plt.axes([0.12,0.55,0.35,0.35]))
    ax.append(plt.axes([0.62,0.55,0.35,0.35]))
    ax.append(plt.axes([0.12,0.10,0.35,0.35]))
    ax.append(plt.axes([0.62,0.10,0.35,0.35]))
    
    for iH in range(heights.size):
    #for iH in range(1):
        
        # isolate parameters for that height
        ht = heights[iH]
        idx_ht = np.where(screen[:,htCol]==ht)[0]
        parms_ht = screen[idx_ht,:]
        
        # calculate empirical CDF
        n_recs = parms_ht.shape[0]
        F = np.arange(1,n_recs+1)/(n_recs+1.)
            
        # plot CDFs
        for iP in range(len(parms)):
    
            # data
            x = np.sort(parms_ht[:,fields.index(parms[iP])])
            
            # fit CDF
            F_fit = jr.compositeCDF(x,*p_parms_opt[iP][iH][:-1])
            
            # plot
            ax[iP].plot(x,F_fit-F,label='{:.0f} m'.format(ht))
            ax[iP].set_title(parms[iP])
            ax[iP].set_ylim([-0.05,0.05])
    
        
    # make things pretty
    ax[2].set_xscale('log')
    ax[0].legend(loc=4,fontsize='small')
    plt.suptitle('Dataset: {:s} ({:s})'.format(dataset,dist_type),fontsize='large')

    iPlot += 1



