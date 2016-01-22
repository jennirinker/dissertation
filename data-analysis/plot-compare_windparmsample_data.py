"""
Draw correlated samples from fit distributions and compare to data 
distributions
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
import JR_Library.main as jr

import numpy as np
import os, json
import matplotlib.pyplot as plt
import scipy.io as scio
import scipy.stats

# choose which dataset
#datasets = ['NREL','fluela','PM06']
datasets = ['PM06']
dist_type = 'comp'
iH        = 0
ParmSample = ['Mean_Wind_Speed', 'Sigma_u', 'Tau_u', 'Concentration_u']
NumSamps = 10000

# define directory where wind parameters are stored (unused for matlab)
BaseDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                  'processed_data'

# plot style
plt.style.use(jr.stylepath('duke_paper'))
FigNum = 1
FigSize = (6.0,4.0)
colors = ['#235F9C', '#C0504D', '#F79646', '#8064A2', \
            '#9BBB59', '#8C564B', '#17BECF', '#BCBD22', '#7F7F7F', '#E377C2', \
            '#262626', '#FFD960']
            
# -----------------------------------------------------------------------------

F_samp = (1. + np.arange(NumSamps))/(1. + NumSamps)

iPlot = 0
for dataset in datasets:

    # ---------- load data ----------
    
    print('\nPlotting {:s} marginal '.format(dist_type) + \
            'distributions to dataset \"{:s}\"\n'.format(dataset))
    
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
    heights = jr.datasetSpecs(dataset_flag)['IDs']
    htCol  = fields.index('ID')
    
    WindParms = jr.SampleWindParameters(NumSamps,dataset,BaseDir,ParmSample,iH,
                         URange=[-float('inf'),11])
    
    # ---------- plot distributions ----------
    
    # initialize figure
    plt.figure(FigNum+iPlot,figsize=FigSize)
    plt.clf()
    ax = []
    ax.append(plt.axes([0.12,0.55,0.35,0.35]))
    ax.append(plt.axes([0.62,0.55,0.35,0.35]))
    ax.append(plt.axes([0.12,0.10,0.35,0.35]))
    ax.append(plt.axes([0.62,0.10,0.35,0.35]))
    
    # isolate parameters for that height
    ht = heights[iH]
    idx_ht = np.where(clean[:,htCol]==ht)[0]
    parms_ht = clean[idx_ht,:]
    
    # calculate empirical CDF
    n_recs = parms_ht.shape[0]
    F_dat = np.arange(1,n_recs+1)/(n_recs+1.)
        
    # plot CDFs
    for iP in range(len(ParmSample)):

        # data
        x = np.sort(parms_ht[:,fields.index(ParmSample[iP])])
        
        # plot
        ax[iP].plot(x,F_dat,label='Empirical'.format(ht))
        ax[iP].plot(np.sort(WindParms[:,iP]),F_samp,
                    label='Sampled'.format(ht))
        ax[iP].set_ylim([0,1])
        
        poop
    
        
    # make things pretty
    ax[2].set_xscale('log')
    ax[0].legend(loc=4,fontsize='small')
    plt.suptitle('Dataset: {:s} ({:s})'.format(dataset,dist_type),
                 fontsize='large')

    iPlot += 1
