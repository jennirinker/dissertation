"""
Compare the residuals before and after generalized pareto distributions
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
import JR_Library.main as jr

import numpy as np
import os, json
import matplotlib.pyplot as plt
import string

# =========================== user variables ==================================

# choose which dataset
#dataset,fignum = 'NREL-mat', 1
#dataset,fignum = 'NREL', 2
#dataset,fignum = 'fluela', 3
dataset,fignum = 'PM06', 4

# define directory where wind parameters are stored (unused for matlab)
basedir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                  'processed_data'

# plot style
plt.style.use(jr.stylepath('duke_paper'))
FigSize = (6.0,3.5)

colors = ['#235F9C', '#C0504D', '#F79646', '#8064A2', \
            '#9BBB59', '#8C564B', '#17BECF', '#BCBD22', '#7F7F7F', '#E377C2', \
            '#262626', '#FFD960']

dx, dy = 0.50, 0.48
xPlot = np.arange(0.12,0.9,dx)
yPlot = np.arange(0.15,0.9,dy)[::-1]
w,h = dx*0.70, dy*0.57

xlabels = ['Mean Wind Speed (m/s)','Turbulence (m/s)',\
           'Kaimal Time Scale (s)','Concentration Parameter']

# ============================ load data ======================================

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

# ========================== plot residuals ===============================

# initialize figure
plt.figure(fignum,figsize=FigSize)
plt.clf()
ax = []
ax.append(plt.axes([xPlot[0],yPlot[0],w,h]))
ax.append(plt.axes([xPlot[1],yPlot[0],w,h]))
ax.append(plt.axes([xPlot[0],yPlot[1],w,h]))
ax.append(plt.axes([xPlot[1],yPlot[1],w,h]))

labels = string.ascii_letters
for idist in range(2):
    dist_type = ['sing','comp'][idist]

    # load distribution information
    dist_fname = '{:s}_6dist_{:s}_parms.txt'.format(dataset,dist_type)
    dist_fpath = os.path.join(basedir,dist_fname)
    with open(dist_fpath,'r') as f:
        dist_dict  = json.load(f)
    p_parms_opt = dist_dict['p_parms_opt']
    parms       = dist_dict['parms']
    
    parms[2] = 'Tau_u'
    
    for iH in range(heights.size):
    #for iH in range(1):
        
        # isolate parameters for that height
        ht = heights[iH]
        idx_ht = np.where(clean[:,htCol]==ht)[0]
        parms_ht = clean[idx_ht,:]
        
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
            ax[iP].plot(x,F_fit-F,color=colors[idist],
                        label='{:.0f} m'.format(ht))
                        
# make things pretty
ax[2].set_xscale('log')
for iP in range(len(parms)):
    ax[iP].set_xlabel(xlabels[iP],
                      fontsize='small')
    ax[iP].set_ylabel('CDF Residual',
                      fontsize='small')
    ylim = np.abs(ax[iP].get_ylim()).max()
    ax[iP].set_ylim([-ylim,ylim])
    ax[iP].locator_params(axis='y',nbins=4)
    ax[iP].text(0.88,0.82,'({:s})'.format(labels[iP]),
                fontsize='medium',
                transform=ax[iP].transAxes)




