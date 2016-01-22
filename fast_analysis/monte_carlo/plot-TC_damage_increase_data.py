"""
Calculate the decrease in lifetime for IEC caused by temporal coherence
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
import os, pickle, string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# plot style
plt.style.use(jr.stylepath('duke_paper'))

# wind datasets
datasets = ['NREL','fluela','PM06']
DSNames = ['M4','Fluela','Plaine Morte','Texas Tech']

# define turbine name and run name
#TurbNames = ['WP0.75A08V00','WP1.5A08V03',
#              'WP3.0A02V02','WP5.0A04V00']
#RunNames = ['Peregrine','TestRun','SmallRun','BigRun2']
TurbName = 'WP5.0A04V00'
RunName  = 'BigRun2'

parameters = [['DEL-h','RootMFlp1','MN-m',1000.,12],
              ['DEL-h','HSShftTq','kN-m',1,5],
              ['DEL-h','TwrBsMyt','MN-m',1000,5]]

FigNum = 1
FigSize = (6.0,3.5)
Labels  = ['$\\rho=0$','$\\rho\\neq0$']

savefig = 0

# base directory where the stats are stored
BaseDmgDictDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                'processed_data'
SaveDir = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
            '2016-02-15_dissertation\\figures'

colors = ['#235F9C', '#C0504D', '#F79646', '#8064A2', \
            '#9BBB59', '#8C564B', '#17BECF']
markers = ['x','o','^','d','+']
labels = string.ascii_letters

dx1, dy1 = 0.99, 0.96/3.
xPlot1 = np.arange(0.09,0.9,dx1)
yPlot1 = np.arange(0.14,0.9,dy1)[::-1]
w1, h1  = dx1*0.75, dy1*0.60

dx2, dy2 = 0.92, 0.96/3.
xPlot2 = np.arange(0.14,0.9,dx2)
yPlot2 = np.arange(0.10,0.9,dy2)[::-1]
w2, h2  = dx2*0.75, dy2*0.50

# -----------------------------------------------------------------------------

iPlot = 0
DPercIncrs = np.empty((len(parameters),len(datasets)))
TDecrs     = np.empty((len(parameters),len(datasets)))
for i_ds in range(len(datasets)):
    dataset = datasets[i_ds]
    
    # load damage dictionary
    DmgDictName = 'DmgIncr_{:s}.txt'.format(dataset)
    DmgDictPath = os.path.join(BaseDmgDictDir,DmgDictName)
    with open(DmgDictPath,'rb') as DictFile:
        OutDict = pickle.load(DictFile)
    
    # add dictionary keys to local variables
    locals().update(OutDict)
    
    # calculate and save percent decreases
    DPercIncrs[:,i_ds] = (100. * (Fs[:,1] - Fs[:,0]) / Fs[:,0])
    TDecrs[:,i_ds]     = 100. - 100. * (1. + DPercIncrs[:,i_ds]/100.)**(-1)
    
    # ----------------- PLOT 1 -- load distributions ------------------------------
    
    # initialize figure, define figure labels
    fig1 = plt.figure(FigNum+iPlot,figsize=FigSize)
    plt.clf()
    
    majorFormatter = FormatStrFormatter('%.2f')
    for istat in range(len(parameters)):
        stat,parm,units,scale,m = parameters[istat]
    
        # initialize axes
        ax = plt.axes([xPlot1[0],yPlot1[istat],w1,h1])
        
        # loop through rho values
        ymax = 0
        for iRho in range(2):
            
            n = np.concatenate((ns[istat,iRho],[0]))
            b = bs[istat,iRho]
            plt.step(b,n,where='post',
                     label=Labels[iRho])
            ymax = max(ymax,n.max()*1.25)
            
        plt.xlabel('{:s} {:s} ({:s})'.format(stat,parm,units),fontsize='small')
        ax.set_ylabel('Frequency',fontsize='small')
        ax.set_ylim([0,ymax])
        plt.locator_params(nbins=5)
        ax.text(0.03,0.75,'({:s})'.format(labels[istat]),transform=ax.transAxes)
        ax.yaxis.set_major_formatter(majorFormatter)
    
        if istat == 0:
            plt.legend(bbox_to_anchor=(1.05, 1),
                       loc=2,borderaxespad=0.,
                       fontsize='x-small')
    
    if savefig:
        FigTitle = 'plot-TC_damage_increase_IEC_hist.eps'
        SavePath = os.path.join(SaveDir,FigTitle)
        fig1.savefig(SavePath)
        print('Figure {:s} saved.'.format(FigTitle))
        
    iPlot += 1
    
    
# plot percent 
fig2 = plt.figure(FigNum+iPlot,figsize=(FigSize))
plt.clf()

for istat in range(len(parameters)):
    ax = plt.axes([xPlot2[0],yPlot2[istat],w2,h2])
    
    for i_ds in range(len(datasets)):
        
        perc_incr = DPercIncrs[istat,i_ds]
        plt.scatter(i_ds,perc_incr,lw=1,
                    marker=markers[0],facecolors='none',
                    edgecolors=colors[0])
                       
#        ax.set_ylim([-10,np.ceil(1.2*DPercIncrs.max()*10)/10.])
#        ax.set_xlim([-1.5,iVref*(2+len(Irefs)) + iIref+1])
    ax.set_xticks(range(len(datasets)))
    ax.set_xticklabels(DSNames[:len(datasets)],
                        fontsize='x-small')
    ax.text(0.02,0.70,'({:s})'.format(labels[istat]),transform=ax.transAxes)
    plt.locator_params(axis='y',nbins=5)
    ax.yaxis.grid()
    ax.tick_params('both', length=0, which='major')    
    ax.set_title('{:s}'.format(parameters[istat][1]),fontsize='small')
    ax.set_ylabel('% Increase\nLifetime\nDamage',fontsize='small')
    
        
if savefig:
    FigTitle = 'plot-TC_damage_increase_IEC_perc.eps'
    SavePath = os.path.join(SaveDir,FigTitle)
    fig2.savefig(SavePath)
    print('Figure {:s} saved.'.format(FigTitle))
        
    iPlot += 1