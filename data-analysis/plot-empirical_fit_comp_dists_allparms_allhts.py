"""
Plot empirical and fit distributions and errors for all parameters and heights
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
import JR_Library.main as jr

import scipy.io as scio
import numpy as np
import json, os
import matplotlib.pyplot as plt
import matplotlib as mpl

# save figure?
saveFig = 0

# choose which dataset
#dataset = 'NREL-mat'
dataset = 'NREL'
#dataset = 'fluela'
#dataset = 'PM06'

# define directory where wind parameters are stored (unused for matlab)
basedir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                  'processed_data'

# plot parameters
plt.style.use(jr.stylepath('duke_paper'))
colors = ['#235F9C', '#C0504D', '#F79646', '#8064A2', '#9BBB59', \
    '#8C564B', '#17BECF', '#BCBD22', '#7F7F7F', '#E377C2']
xPos   = [0.04,0.46]
yPos1  = [0.56,0.05]
yPos2  = [0.76,0.25]
    
# ======================= load data/parameters ==============================

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

# load fit parameters
fname = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
        'processed_data\\NREL_optdist_comp_parms.txt'
with open(fname,'r') as f:
    p_parms_opt = json.load(f)

# ============================= plot data ===================================

# initialize figure
fig = plt.figure(1,figsize=(15,7))
plt.clf()

# empirical CDF
N     = clean.shape[0]
F_emp = np.arange(N)/(N+1.)

# loop through parameters
for iP in range(4):
    
    # get axes indices, overwrite default
    iX = (iP % 2)
    iY = int(iP/2.)
    ax2 = plt.axes([xPos[iX],yPos1[iY],0.35,0.18])
    ax1 = plt.axes([xPos[iX],yPos2[iY],0.35,0.18])
    ax1.xaxis.set_tick_params(length=4)
    ax2.xaxis.set_tick_params(length=4)
    ax1.yaxis.set_tick_params(length=4)
    ax2.yaxis.set_tick_params(length=4)
    
    # plot dotted lines to mark zero error
    if iP == 0:
        ax2.plot([0,30],[0,0],'k:',lw=2)
        ax2.set_xticklabels(['0','5','10','15','20','25','30'])
    elif iP == 1: 
        ax2.plot([0,7],[0,0],'k:')
        ax2.set_xticklabels(['0','1','2','3','4','5','6','7'])
    elif iP == 2:
        ax2.plot([10**(-1),10**(5)],[0,0],'k:',lw=2)
        ax1.set_xscale('log')
        ax2.set_xscale('log')
    elif iP == 3:
        ax1.plot([0.027,0.027],[0,1],'k--',lw=2)
        ax2.plot([0,1],[0,0],'k:',lw=2)
        ax2.set_xticklabels(['0.0','0.2','0.2','0.6','0.8','1.0'])
    
    # loop through heights
    for iH in range(6):
                
        # extract data
        x = clean[:,5*iH + iP]
        x = np.sort(x)
        
        # calculate fit CDFs
        p_comp = p_parms_opt[iP][iH]
        F = jr.compositeCDF(x,p_comp[0],p_comp[1],p_comp[2],p_comp[3])
        
        # Lower plot: CDF error
        CDF_err = F - F_emp
        ax2.plot(x,CDF_err,color=colors[iH],lw=2)
        
        # upper plot: CDF
        ax1.plot(x,F_emp,color=colors[iH],lw=2)
        skip = 800
        ax1.plot(x[::skip],F[::skip],'^',mec=colors[iH],mfc='None',mew=2)
        
    if (iP == 1):
        legStr = ['15 m - Data','15 m - Model','30 m - Data','30 m - Model', \
            '50 m - Data','50 m - Model','76 m - Data','76 m - Model', \
            '100 m - Data','100 m - Model','131 m - Data','131 m - Model']
        plt.legend(legStr, fontsize = 'large', \
            bbox_to_anchor=(1.12, 1), loc=2, borderaxespad=0.)        

        
    # make axes look good
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_xticklabels([''])
    ax1.set_yticklabels(['','0.2','0.4','0.6','0.8','1.0'])
    
    
    ax2.set_ylim([-0.02,0.02])
    ax2.set_yticks([-0.02,-0.015,-0.01,-0.005,0,0.005,0.01,0.015,0.02])
    ax2.set_yticklabels(['-0.02','','-0.01','','0.00','','0.01','','0.02'])

    
    ax1.text(0.97,0.10,'CDF', \
            transform=ax1.transAxes,ha='right',va='bottom')
    ax2.text(0.97,0.90,'CDF Error', \
            transform=ax2.transAxes,ha='right',va='top')

plt.text(0.215,0.97,'Mean Wind Speed $U$',fontsize='x-large', \
    transform=fig.transFigure,ha='center',va='center')
plt.text(0.635,0.97,'10-Minute Standard Deviation $\sigma_u$',fontsize='x-large', \
    transform=fig.transFigure,ha='center',va='center')
plt.text(0.215,0.46,'Kaimal Length Scale $L$',fontsize='x-large', \
    transform=fig.transFigure,ha='center',va='center')
plt.text(0.635,0.46,r'Concentration Parameter $\rho$',fontsize='x-large', \
    transform=fig.transFigure,ha='center',va='center')   
   
# ADD TABLE TO PLOT ('usetex' flag messes up settings)
#mpl.rc('text',usetex='True')
#tabStr = r'\begin{tabular}{ c | c c c c} & $U$ & $\sigma_u$ & $\rho$ & $L$ \\\hline' \
#    + r' 15 m & C2 & GG & GEV & GG \\ 30 m & C2 & GG & GEV & B \\' + \
#    r'50 m & EW & EW & GEV & GG \\ 76 m & C2 & GG & GEV & LN \\' + \
#    r'100 m & C2 & GG & GEV & GG \\ 131 m & LN & EW & GEV & GG\end{tabular}'
#plt.text(0.91,0.42,'Distribution models:',fontsize='medium', \
#    transform=fig.transFigure,ha='center',va='bottom')  
#plt.text(0.91,0.40,tabStr,fontsize='medium', \
#    transform=fig.transFigure,ha='center',va='top')  
#defStr = r'\begin{tabular}{rl} C2:&Chi-Squared\\GG:&Generalized Gamma\\' + \
#    r'GEV:&Generalized Extreme Value\\B:&Beta\\EW:&Exponential Weibull\\'+ \
#    r'LN:&Lognormal\end{tabular}'
#plt.text(0.85,0.17,defStr,fontsize='small', \
#    transform=fig.transFigure,ha='left',va='top')   

# save figure if requested
if saveFig:
    plt.savefig('C:\\Users\\jrinker\\Dropbox\\presentations\\' + \
        '2015-06-10_nawea\\presentation\\pictures\\15-all_marginals.png')
    print 'Figure saved.'

