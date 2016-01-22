""" Plot correlation plots
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import numpy as np
import scipy.stats
import JR_Library.main as jr
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.io as scio
import os,string

# datasets to plot correlations for
datasets = ['NREL','fluela','PM06']
dataset = 'NREL'

# plot options
plt.style.use(jr.stylepath('duke_paper'))
SaveFig = 0
SaveDir = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
            '2016-02-15_dissertation\\figures'
p_iHplot = 1          # height to plot
FigNum = 1

# axes locations
xPos  = [0.09,0.33,0.57,0.81]
yPos  = xPos[::-1]
w     = 0.15
h     = 0.15

parm_names = ['$U$','$\\sigma_u$','$\\tau_u$','$\\rho_u$','$\\zeta$']
parms = ['Mean_Wind_Speed','Sigma_u','Tau_u','Concentration_u','MO_Length_virt']

# define directory where wind parameters are stored (unused for matlab)
BaseDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                  'processed_data'

# -----------------------------------------------------------------------------

iPlot = 0
for dataset in datasets:

    # load file with all heights for plotting histograms
    AllHtsName = '{:s}-metadata_allheights.mat'.format(dataset)
    AllHtsPath = os.path.join(BaseDir,AllHtsName)
    outdict = scio.loadmat(AllHtsPath,squeeze_me=True)
    values  = outdict['values']
    IDs = outdict['IDs']
    all_fields  = [s.rstrip() for s in outdict['all_fields']]
    
    iHplot = (len(IDs)-1)*p_iHplot
    
    # load correlations
    DictName   = '{:s}_correlations.mat'.format(dataset)
    DictPath   = os.path.join(BaseDir,DictName)
    CorrDict   = scio.loadmat(DictPath,squeeze_me=True)
    IDs,R      = CorrDict['IDs'], CorrDict['R']
    all_fields = [s.rstrip() for s in CorrDict['all_fields']]
    
    # initialize figure
    fig = plt.figure(FigNum+iPlot,figsize=(4.5,4.5))
    plt.clf()
    
    
    labels = string.ascii_letters
    
    # standardize the random variables
    F = (np.arange(values.shape[0])+1.)/(values.shape[0]+1.)
    gvals = np.empty(values.shape)[:,1:]
    for i_parm in range(gvals.shape[1]):
        gvals[np.argsort(values[:,i_parm+1]),i_parm] = scipy.stats.norm.ppf(F)
    
    
    # loop through first variable
    iAx = 0
    for iParm1 in range(len(parms)):
        parm1 = parms[iParm1]
        i1    = iHplot*len(all_fields) + all_fields.index(parm1)
        G_x   = gvals[:,i1]
        
        # loop through second variable
        for iParm2 in range(iParm1+1,len(parms)):
            parm2 = parms[iParm2]
            i2    = iHplot*len(all_fields) + all_fields.index(parm2)
            G_y   = gvals[:,i2]
            
            # calculate the pearson coefficient
            r = R[i1,i2]
            
            # make axes
            ax = plt.axes([xPos[iParm2-1],yPos[iParm1],w,h])  
              
            # plot
            plt.hist2d(G_x,G_y,bins=30,cmap=plt.get_cmap('Blues'))
            
            # make axes pretty
            plt.axis('equal')
            jr.removeSpines(ax)
            ax.axes.get_xaxis().set_ticks([-4,-3,-2,-1,0,1,2,3,4])
            ax.axes.get_xaxis().set_ticklabels(['-4','','-2','','0','','2','','4'],
                                               fontsize='x-small')
            ax.axes.get_yaxis().set_ticks([-4,-3,-2,-1,0,1,2,3,4])
            ax.axes.get_yaxis().set_ticklabels(['-4','','-2','','0','','2','','4'],
                                               fontsize='x-small')
            
            # add plot annotations
            ax.text(0.23,0.88,'({:s})'.format(labels[iAx]),ha='center', \
                    va='center',fontsize='small',transform=ax.transAxes)
            ax.text(0.98,0.15,'r = {:.2f}'.format(r),ha='right', \
                    va='center',fontsize='small',transform=ax.transAxes)
            ax.text(-6.3,0,'{}'.format(parm_names[iParm1]),ha='center', rotation='vertical', \
                    va='center',fontsize='medium')
            ax.text(0,-7.1,'{}'.format(parm_names[iParm2]),ha='center', \
                    va='center',fontsize='medium')
    
            iAx += 1
        
    if SaveFig:
        FigTitle = 'plot-{:s}_correlations.eps'.format(dataset)
        SavePath = os.path.join(SaveDir,FigTitle)
        fig.savefig(SavePath)
        print('Figure {:s} saved.'.format(FigTitle))
        
    iPlot += 1
    
    