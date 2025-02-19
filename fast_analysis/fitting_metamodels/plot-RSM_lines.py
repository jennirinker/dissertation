"""
Fitting polynomial surface to FAST statistics
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
import os, pickle
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# plot style
plt.style.use(jr.stylepath('duke_paper'))

# define turbine name and run name
#TurbNames = ['WP0.75A08V00','WP1.5A08V03',
#              'WP3.0A02V02','WP5.0A04V00']
#RunNames = ['Peregrine','TestRun','SmallRun','BigRun2']
TurbNames = ['WP5.0A04V00']
RunName  = 'BigRun2'

FigNum = 4
FigSize = (6.0,7.5)

#parameters = [['max','RootMFlp1','MN-m',1000.]]
parameters = [['max','RootMFlp1','MN-m',1000.],
              ['DEL-h','RootMFlp1','MN-m',1000.],
              ['max','HSShftTq','kN-m',1],
              ['DEL-h','HSShftTq','kN-m',1],
              ['max','TwrBsMyt','MN-m',1000],
              ['DEL-h','TwrBsMyt','MN-m',1000],
              ['mean','GenPwr','MW',1000.]]

# base directory where the stats are stored
BaseStatDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                'processed_data\\proc_stats'

# -----------------------------------------------------------------------------

# get wind parameters for that run
WindParms = jr.RunName2WindParms(RunName)
URefs, Is, Ls, rhos, n_dups = WindParms['URefs'],WindParms['Is'], \
                              WindParms['Ls'],WindParms['rhos'], \
                              WindParms['n_dups']
WindParms = [URefs,Is,np.log10(Ls),rhos]
UPlot     = np.linspace(3,22,201)

# loop through turbines and stats
iPlot = 0
for TurbName in TurbNames:
    
    TurbRSMDictName = '{:s}_RSM.bdat'.format(TurbName)
    TurbRSMDictPath = os.path.join('RSMs',TurbRSMDictName)
    with open(TurbRSMDictPath,'rb') as f:
        TurbRSMDict = pickle.load(f)
    
    for stat,parm,units,scale in parameters:
        
        # load the stats data
        x, y = jr.LoadFASTStats(RunName,TurbName,stat,parm,
                                scale=scale)
        
        # load the RSM data
        DictKey = '{:s}_{:s}'.format(parm,stat)
        RSMDict = TurbRSMDict[DictKey]
        ps, cs  = RSMDict['ps_red'], RSMDict['cs_red']
        
        # ================= plot data vs I =====================

        # initialize figure
        fig = plt.figure(FigNum+iPlot,figsize=FigSize)
        plt.clf()
                
        # plot data
#        ymax = max(-np.floor(yplot.min()), np.ceil(yplot.max()))
        for iL in range(len(WindParms[2])):
            logL = WindParms[2][iL]
            for iRho in range(len(WindParms[3])):
                rho = WindParms[3][iRho]
                
                # create axes
                iplot = iRho*len(WindParms[2]) + iL + 1
                ax = fig.add_subplot(len(WindParms[3]),len(WindParms[2]),iplot)
                
                # plot data
                for I in np.arange(0,0.6,0.1):
                    x = np.empty((len(UPlot),4))
                    x[:,0],x[:,1],x[:,2],x[:,3] = UPlot,I,logL,rho
                    Xv = jr.myvander(x,ps)
                    yplot = np.dot(Xv,cs)
                    
                    ax.plot(UPlot,yplot,
                            label='$I$ = {:.1f}'.format(WindParms[1][iI]))
                   
                # prettify axes
#                ax.set_ylim([-ymax,ymax])
                ax.set_xlim([4,24])
                plt.locator_params(axis='x',nbins=5)
                jr.removeSpines(ax)
                
                # put x and y labels on left column and bottom row only
                if iRho < len(WindParms[3])-1:
                    ax.set_xticklabels([])
                if iL > 0:
                    ax.set_yticklabels([])
                    
                
                # create legend
                if (iL == 0) and (iRho == 0):
                    plt.legend(bbox_to_anchor=(-0.10, 0.94, 4.0, 0.94),
                               loc=3,ncol=5)
        
        # scale subplots and add text labels
        xbord,ybord = 0.07,0.025
        plt.tight_layout(rect=[xbord,ybord,1.01,0.94])
        plt.figtext(0.5+xbord/2.,0.02,'Mean Wind Speed [m/s]',
                    ha='center',va='center')
        plt.figtext(0.07,0.5+ybord/2.,'{:s} {:s} [{:s}]'.format(stat,parm,units),
                    ha='center',va='center',rotation='vertical')
        
        for iRho in range(len(WindParms[3])):
            plt.figtext(0.025,0.85-0.88*iRho/len(WindParms[3]),
                        r'$\rho$ = {:.1f}'.format(WindParms[3][iRho]),
                        ha='center',va='center',rotation='vertical')
        for iL in range(len(WindParms[2])):
            plt.figtext(0.23 + 0.88*iL/len(WindParms[2]),0.98,
                        r'log$_{10}$L = ' + '{:.1f}'.format(WindParms[2][iL]),
                        ha='center',va='center')
                        
        iPlot += 1
        
        
        poop1
        # ================= plot data vs rho =====================
        
        # initialize figure
        fig = plt.figure(FigNum+iPlot,figsize=FigSize)
        plt.clf()
        
        # plot data
#        ymax = max(-np.floor(yplot.min()), np.ceil(yplot.max()))
        for iL in range(len(WindParms[2])):
            logL = WindParms[2][iL]
            for iI in range(len(WindParms[1])):
                I = WindParms[1][iI]
                
                # create axes
                iplot = iI*len(WindParms[2]) + iL + 1
                ax = fig.add_subplot(len(WindParms[1]),len(WindParms[2]),iplot)
                
                # plot data
                for iRho in range(len(WindParms[3])):
                    rho = WindParms[3][iRho]
                    x = np.empty((len(UPlot),4))
                    x[:,0],x[:,1],x[:,2],x[:,3] = UPlot,I,logL,rho
                    Xv = jr.myvander(x,ps)
                    yplot = np.dot(Xv,cs)
                    
                    ax.plot(UPlot,yplot,
                            label=r'$\rho$ = {:.1f}'.format(WindParms[3][iRho]))
                    
                # prettify axes
#                ax.set_ylim([-ymax,ymax])
                ax.set_xlim([4,24])
                plt.locator_params(axis='x',nbins=5)
                jr.removeSpines(ax)
                
                # put x and y labels on left column and bottom row only
                if iRho < len(WindParms[1])-1:
                    ax.set_xticklabels([])
                if iL > 0:
                    ax.set_yticklabels([])
                
                # create legend
                if (iL == 0) and (iI == 0):
                    plt.legend(bbox_to_anchor=(-0.10, 0.94, 4.0, 0.94),
                               loc=3,ncol=5)
        
        # scale subplots and add text labels
        xbord,ybord = 0.07,0.025
        plt.tight_layout(rect=[xbord,ybord,1.01,0.94])
        plt.figtext(0.5+xbord/2.,0.02,'Mean Wind Speed [m/s]',
                    ha='center',va='center')
        plt.figtext(0.07,0.5+ybord/2.,'{:s} {:s} [{:s}]'.format(stat,parm,units),
                    ha='center',va='center',rotation='vertical')
        
        for iI in range(len(WindParms[1])):
            plt.figtext(0.025,0.85-0.88*iI/len(WindParms[1]),
                        r'$I$ = {:.1f}'.format(WindParms[1][iI]),
                        ha='center',va='center',rotation='vertical')
        for iL in range(len(WindParms[2])):
            plt.figtext(0.23 + 0.88*iL/len(WindParms[2]),0.98,
                        r'log$_{10}$L = ' + '{:.1f}'.format(WindParms[2][iL]),
                        ha='center',va='center')
                        
        iPlot += 1

        # ================= plot data vs L =====================
        
        # initialize figure
        fig = plt.figure(FigNum+iPlot,figsize=FigSize)
        plt.clf()
        
        # plot data
#        ymax = max(-np.floor(yplot.min()), np.ceil(yplot.max()))
        for iI in range(len(WindParms[1])):
            I = WindParms[1][iI]
            for iRho in range(len(WindParms[3])):
                rho = WindParms[3][iRho]
                
                # create axes
                iplot = iI*len(WindParms[3]) + iRho + 1
                ax = fig.add_subplot(len(WindParms[1]),len(WindParms[3]),iplot)
                
                # plot data
                for iL in range(len(WindParms[2])):
                    logL = WindParms[2][iL]
                    x = np.empty((len(UPlot),4))
                    x[:,0],x[:,1],x[:,2],x[:,3] = UPlot,I,logL,rho
                    Xv = jr.myvander(x,ps)
                    yplot = np.dot(Xv,cs)
                    
                    ax.plot(UPlot,yplot,
                            label='$L$ = $10^{' + '{:.1f}'.format(WindParms[2][iL]) + \
                                '}$')
                    
                # prettify axes
#                ax.set_ylim([-ymax,ymax])
                ax.set_xlim([4,24])
                plt.locator_params(axis='x',nbins=5)
                jr.removeSpines(ax)
                
                # put x and y labels on left column and bottom row only
                if iI < len(WindParms[1])-1:
                    ax.set_xticklabels([])
                if iRho > 0:
                    ax.set_yticklabels([])
                
                # create legend
                if (iI == 0) and (iRho == 0):
                    plt.legend(bbox_to_anchor=(-0.10, 0.94, 4.0, 0.94),
                               loc=3,ncol=4)
        
        # scale subplots and add text labels
        xbord,ybord = 0.07,0.025
        plt.tight_layout(rect=[xbord,ybord,1.01,0.94])
        plt.figtext(0.5+xbord/2.,0.02,'Mean Wind Speed [m/s]',
                    ha='center',va='center')
        plt.figtext(0.07,0.5+ybord/2.,'{:s} {:s} [{:s}]'.format(stat,parm,units),
                    ha='center',va='center',rotation='vertical')
        
        for iI in range(len(WindParms[1])):
            plt.figtext(0.025,0.85-0.88*iI/len(WindParms[1]),
                        r'$I$ = {:.1f}'.format(WindParms[1][iI]),
                        ha='center',va='center',rotation='vertical')
        for iL in range(len(WindParms[2])):
            plt.figtext(0.23 + 0.88*iL/len(WindParms[2]),0.98,
                        '$\\rho$ = ' + '{:.1f}'.format(WindParms[3][iRho]),
                        ha='center',va='center')
                        
        iPlot += 1
        
        poop

