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

# define turbine name and run name
#TurbNames = ['WP0.75A08V00','WP1.5A08V03',
#              'WP3.0A02V02','WP5.0A04V00']
#RunNames = ['Peregrine','TestRun','SmallRun','BigRun2']
TurbName = 'WP5.0A04V00'
RunName  = 'BigRun2'

FigNum = 1
FigSize = (6.0,3.5)

parameters = [['max','RootMFlp1','MN-m',1000.],
              ['max','HSShftTq','kN-m',1],
              ['max','TwrBsMyt','MN-m',1000]]

DmgDictName = 'UltIncr_IEC.txt'

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

dx1, dy1 = 0.96, 0.96/3.
xPlot1   = np.arange(0.11,0.9,dx1)
yPlot1   = np.arange(0.14,0.9,dy1)[::-1]
w1, h1   = dx1*0.75, dy1*0.60

dx2, dy2 = 0.92, 0.96/3.
xPlot2   = np.arange(0.14,0.9,dx2)
yPlot2   = np.arange(0.10,0.9,dy2)[::-1]
w2, h2   = dx2*0.75, dy2*0.50

# -----------------------------------------------------------------------------

# load damage dictionary
DmgDictPath = os.path.join(BaseDmgDictDir,DmgDictName)
with open(DmgDictPath,'rb') as DictFile:
    OutDict = pickle.load(DictFile)

# add dictionary keys to local variables
locals().update(OutDict)

# ----------------- PLOT 1 -- load distributions ------------------------------

# WT class 1, turbulence class A
iVref, iIref = 0, 0

# initialize figure, define figure labels
fig1 = plt.figure(FigNum,figsize=FigSize)
plt.clf()

# loop through statistics
iPlot = 0

majorFormatter = FormatStrFormatter('%.2f')
for istat in range(len(parameters)):
    stat,parm,units,scale = parameters[istat]

    # initialize axes
    ax = plt.axes([xPlot1[0],yPlot1[istat],w1,h1])
    
    # loop through rho values
    ymax = 0
    for iRho in range(len(rhos)):
        rho = rhos[iRho]
        
        n = np.concatenate((ns[istat,iRho,iVref,iIref],[0]))
        b = bs[istat,iRho,iVref,iIref]
        plt.step(b,n,where='post',
                 label='$\\rho={:.1f}$'.format(rho))
        ymax = max(ymax,n.max()*1.25)
        
    # plot sampled distribution manually
    n = np.concatenate((ns[istat,len(rhos),iVref,iIref],[0]))
    b = bs[istat,len(rhos),iVref,iIref]
    plt.step(b,n,where='post',lw=2,
             label='M4 $\\rho$'.format(rho))
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
    FigTitle = 'plot-TC_ultimate_increase_IEC_hist.eps'
    SavePath = os.path.join(SaveDir,FigTitle)
    fig1.savefig(SavePath)
    print('Figure {:s} saved.'.format(FigTitle))

# ----------------- PLOT 2 -- percent increases ------------------------------

# calculate percent increases
DPercIncrs = np.empty((len(parameters),len(rhos),len(Vrefs),len(Irefs)))
TDecrs     = np.empty((len(parameters),len(rhos),len(Vrefs),len(Irefs)))
for iVref in range(len(Vrefs)):
    for iIref in range(len(Irefs)):
        
        Dhat0 = np.dot(Fs[:,0,iVref,iIref].reshape(len(parameters),1),
                       np.ones((1,len(rhos)+1)))
        DPercIncr = (100. * (Fs[:,:,iVref,iIref] - Dhat0) / Dhat0)[:,1:]
        TDecr = 100. - 100. * (1. + DPercIncr/100.)**(-1)
        
        DPercIncrs[:,:,iVref,iIref] = DPercIncr
        TDecrs[:,:,iVref,iIref]     = TDecr


# plot percent 
fig2 = plt.figure(FigNum+1,figsize=(FigSize))
plt.clf()

for istat in range(len(parameters)):
    ax = plt.axes([xPlot2[0],yPlot2[istat],w2,h2])
    
    for iVref in range(len(Vrefs)):
        for iIref in range(len(Irefs)):
        
            x = iVref*(2+len(Irefs)) + iIref

            for iRho in range(len(rhos)-1):
                
                perc_incr = DPercIncrs[istat,iRho,iVref,iIref]
                plt.scatter(x,perc_incr,lw=1,
                            marker=markers[iRho],facecolors='none',
                            edgecolors=colors[iRho],
                             label='$\\rho={:.1f}$'.format(rhos[iRho+1]) \
                                 if (iVref+iIref == 0) else "")
                
            perc_incr = DPercIncrs[istat,len(rhos)-1,iVref,iIref]
            plt.scatter(x,perc_incr,lw=1,
                        marker='o',facecolors=colors[len(rhos)],
                        edgecolors=colors[len(rhos)],
                         label='M4 $\\rho$' \
                             if (iVref+iIref == 0) else "")
                       
    ax.set_ylim([-2,np.ceil(1.2*DPercIncrs.max()*10)/10.])
    ax.set_xlim([-1.5,iVref*(2+len(Irefs)) + iIref+1])
    ax.set_xticks([0,1,2,5,6,7,10,11,12])
    ax.set_xticklabels(['\nTC-A','Wind Turbine Class 1\nTC-B','\nTC-C',
                        '\nTC-A','Wind Turbine Class 2\nTC-B','\nTC-C',
                        '\nTC-A','Wind Turbine Class 3\nTC-B','\nTC-C'],
                        fontsize='x-small')
    ax.text(0.02,0.70,'({:s})'.format(labels[istat]),transform=ax.transAxes)
    plt.locator_params(axis='y',nbins=4)
    ax.yaxis.grid()
    ax.tick_params('both', length=0, which='major')    
    ax.set_title('{:s}'.format(parameters[istat][1]),fontsize='small')
    ax.set_ylabel('% Increase\nUltimate\nLoad',fontsize='small')
    
    if istat == 0:
        plt.legend(bbox_to_anchor=(1.05, 1),
                   loc=2,borderaxespad=0.,
                   fontsize='x-small')
        
if savefig:
    FigTitle = 'plot-TC_ultimate_increase_IEC_perc.eps'
    SavePath = os.path.join(SaveDir,FigTitle)
    fig2.savefig(SavePath)
    print('Figure {:s} saved.'.format(FigTitle))