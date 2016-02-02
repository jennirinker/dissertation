"""
Compare the lifetime damage for IEC and data
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
FigSize = (6.0,3.75)

savefig = 0

# wind datasets
datasets = ['NREL','fluela','PM06','texastech']
DSNames = ['IEC:\nWTC-I,\nTC-A','IEC:\nWTC-III\nTC-C','M4','Fluela',
           'Plaine\nMorte','Texas\nTech']

# base directory where the stats are stored
BaseDmgDictDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                'processed_data'
SaveDir = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
            '2016-02-15_dissertation\\figures'

colors = ['#235F9C', '#C0504D', '#F79646', '#8064A2', \
            '#9BBB59', '#8C564B', '#17BECF']
markers = ['x','o','^','d','+']
labels = string.ascii_letters
LabelStr = ['$\\rho=0$','$\\rho\\neq0$']

dx1, dy1 = 0.99, 0.96/3.
xPlot1 = np.arange(0.09,0.9,dx1)
yPlot1 = np.arange(0.14,0.9,dy1)[::-1]
w1, h1  = dx1*0.75, dy1*0.60

dx2, dy2 = 0.92, 0.96/3.
xPlot2 = np.arange(0.16,0.9,dx2)
yPlot2 = np.arange(0.12,0.9,dy2)[::-1]
w2, h2  = dx2*0.75, dy2*0.50

# -----------------------------------------------------------------------------

# initialize array of total forces
Fs_all = np.empty((3,2,6))

# load damage dictionary for IEC
DmgDictName = 'DmgIncr_IEC.txt'
DmgDictPath = os.path.join(BaseDmgDictDir,DmgDictName)
with open(DmgDictPath,'rb') as DictFile:
    OutDict = pickle.load(DictFile)

# add dictionary keys to local variables for IEC
locals().update(OutDict)

# assign F values to total array
Fs_all[:,:,0] = Fs[:,[0,-1],0,0]
Fs_all[:,:,1] = Fs[:,[0,-1],-1,-1]

# loop through datasets
for i_ds in range(len(datasets)):
    dataset     = datasets[i_ds]
    DmgDictName = 'DmgIncr_{:s}.txt'.format(dataset)
    DmgDictPath = os.path.join(BaseDmgDictDir,DmgDictName)
    with open(DmgDictPath,'rb') as DictFile:
        OutDict = pickle.load(DictFile)
    
    # add dictionary keys to local variables for IEC
    locals().update(OutDict)
    
    # assign F values to total array
    Fs_all[:,:,i_ds+2] = Fs

# ----------------- scatter response ------------------------------

# plot percent 
fig = plt.figure(FigNum,figsize=(FigSize))
plt.clf()

for istat in range(len(parameters)):
    ax = plt.axes([xPlot2[0],yPlot2[istat],w2,h2])
    
    Fsnorm = Fs_all[istat]/Fs_all[istat].min()
    for i_ds in range(Fs_all.shape[-1]):
        
        for iRho in range(2):
            plt.scatter(i_ds,Fsnorm[iRho,i_ds],lw=1,
                        marker=markers[iRho],facecolors='none',
                        edgecolors=colors[iRho],
                         label=LabelStr[iRho] \
                             if (i_ds == 0) else "")
                    
#    ax.set_ylim([-10,np.ceil(1.2*DPercIncrs.max()*10)/10.])
    ax.set_xlim([-1,Fs_all.shape[2]])
    ax.set_xticks(range(Fs_all.shape[2]))
    ax.set_xticklabels(DSNames,
                        fontsize='x-small')
    ax.text(0.02,0.70,'({:s})'.format(labels[istat]),transform=ax.transAxes)
    plt.locator_params(axis='y',nbins=5)
    ax.yaxis.grid()
    ax.tick_params('both', length=0, which='major')    
    ax.set_title('{:s}'.format(parameters[istat][1]),fontsize='small')
    ax.text(-0.16,0.5,'Normalized\nLifetime\nDamage',
            transform=ax.transAxes,rotation='vertical',
            va='center',ha='center',fontsize='small')
    
    if istat == 0:
        plt.legend(bbox_to_anchor=(1.05, 1),
                   loc=2,borderaxespad=0.,
                   fontsize='x-small')
        
if savefig:
    FigTitle = 'plot-TC_damage_data_IEC.eps'
    SavePath = os.path.join(SaveDir,FigTitle)
    fig.savefig(SavePath)
    print('Figure {:s} saved.'.format(FigTitle))