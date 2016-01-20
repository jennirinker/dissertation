"""
IEC Monte Carlo for 5.0 MW turbine
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
import os, pickle, json
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as scio

# plot style
plt.style.use(jr.stylepath('duke_paper'))

# define turbine name and run name
#TurbNames = ['WP0.75A08V00','WP1.5A08V03',
#              'WP3.0A02V02','WP5.0A04V00']
#RunNames = ['Peregrine','TestRun','SmallRun','BigRun2']
TurbNames = ['WP5.0A04V00']
RunName  = 'BigRun2'

FigNum = 2
#FigSize = (6.0,8.5)    # fulll page
FigSize = (4.0,5.25)    # half page
LegLoc = (0.03, 0.92, 0.96, .05)
LegFS = 'x-small'
ULim = [4,24]
TickBins = 4

savefig = 0

zRef  = 90.
shear = 0.2
Vref  = 50
Iref  = 0.16
Uref_lo,Uref_hi = 5,22

cov = 0.1

NumSamps = 5000

parameters = [['max','RootMFlp1','MN-m',1000.]]
#parameters = [['DEL-h','RootMFlp1','MN-m',1000.]]
#parameters = [['max','RootMFlp1','MN-m',1000.],
#              ['DEL-h','RootMFlp1','MN-m',1000.],
#              ['max','HSShftTq','kN-m',1],
#              ['DEL-h','HSShftTq','kN-m',1],
#              ['max','TwrBsMyt','MN-m',1000],
#              ['DEL-h','TwrBsMyt','MN-m',1000],
#              ['mean','GenPwr','MW',1000.]]

# base directory where the stats are stored
BaseStatDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                'processed_data\\proc_stats'
SaveDir = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
            '2016-02-15_dissertation\\figures'
BaseTurbDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
                'dissertation\\FAST_models\\FAST7'
RSMDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
                'dissertation\\fast_analysis\\fitting_metamodels\\RSMs'

# -----------------------------------------------------------------------------


rands = np.random.rand(NumSamps)
randn = np.random.normal(size=NumSamps)

rho = 0.0

## get wind parameters for that run
#WindParms = jr.RunName2WindParms(RunName)
#URefs, Is, Ls, rhos, n_dups = WindParms['URefs'],WindParms['Is'], \
#                              WindParms['Ls'],WindParms['rhos'], \
#                              WindParms['n_dups']
#WindParms = [URefs,Is,np.log10(Ls),rhos]
#UPlot     = np.linspace(3,22,201)
#
## loop through turbines and stats
#iPlot = 0
for TurbName in TurbNames:
    
    print('Turbine {:s}'.format(TurbName))
    
    # load RSM dictionary
    TurbRSMDictName = '{:s}_RSM.bdat'.format(TurbName)
    TurbRSMDictPath = os.path.join(RSMDir,TurbRSMDictName)
    with open(TurbRSMDictPath,'rb') as f:
        TurbRSMDict = pickle.load(f)
    
    # load look up table and turbine dictionry
    SSPath = os.path.join(BaseTurbDir,TurbName,'steady_state',
                          '{:s}_SS.mat'.format(TurbName))
    SSDict = scio.loadmat(SSPath,squeeze_me=True)
    SSData = SSDict['SS']
    SSFields = [s.rstrip() for s in SSDict['Fields']]
    TurbDictPath = os.path.join(BaseTurbDir,TurbName,'parameters',
                          '{:s}_Dict.dat'.format(TurbName))
    with open(TurbDictPath,'r') as DictFile:
        HubHeight = json.load(DictFile)['HH']
    logL = np.log10(jr.IEC_Lambda1(HubHeight)*8.1)
        
    Uhub_lo = Uref_lo*(HubHeight/zRef)**shear
    Uhub_hi = Uref_hi*(HubHeight/zRef)**shear
    F_lo    = 1 - np.exp(-np.pi*(Uhub_lo/0.4/Vref)**2)
    F_hi    = 1 - np.exp(-np.pi*(Uhub_hi/0.4/Vref)**2)
    
    rands   = (F_hi - F_lo)*rands + F_lo
    Uhubs = 0.4*Vref*np.sqrt(-np.log(1-rands)/np.pi)
    sig_us = Iref*(0.75*Uhubs+5.6)

    Urefs = Uhubs/(HubHeight/zRef)**shear
    Is    = sig_us/Urefs
        
    x = np.empty((NumSamps,4))
    x[:,0],x[:,1],x[:,2],x[:,3] = Urefs,Is,logL,rho
    
    for stat,parm,units,scale in parameters:
        
        print('  {:s} {:s}'.format(stat,parm))
        
        # load the RSM data
        DictKey = '{:s}_{:s}'.format(parm,stat)
        RSMDict = TurbRSMDict[DictKey]
        ps, cs  = RSMDict['ps_red'], RSMDict['cs_red']
        
        # calculate mean loads
        Xv = jr.myvander(x,ps)
        mean_loads = np.dot(Xv,cs)
        
        # add randomness
        loads = mean_loads + mean_loads*cov*randn
#        loads = mean_loads
        
        print(loads.max())
        
        plt.figure(4)
        plt.clf()
        
        plt.subplot(1,3,1)
        plt.hist(Urefs,bins=50)
        
        plt.subplot(1,3,2)
        plt.hist(Is,bins=50)
        
        plt.subplot(1,3,3)
        plt.hist(loads,bins=50)
        
#        
#        # ================= plot data vs I =====================
#        
#        # initialize figure
#        fig = plt.figure(FigNum,figsize=FigSize)
#        plt.clf()
#        
#        # plot data
#        for iRho in range(len(WindParms[3])):
#            rho = WindParms[3][iRho]
#            for iL in range(len(WindParms[2])):
#                logL = WindParms[2][iL]
#                
#                # create axes
#                iplot = iRho*len(WindParms[2]) + iL + 1
#                ax = fig.add_subplot(len(WindParms[3]),len(WindParms[2]),iplot)
#                
#                # plot steady-state if applicable
#                if stat in ['mean','max']:
#                    U_Refs = np.linspace(min(WindParms[0]),max(WindParms[0]),
#                                         20)
#                    U_hubs = U_Refs*(HubHeight/zRef)**shear
#                    y_SS   = np.interp(U_hubs,SSData[:,SSFields.index('WindVxi')],
#                                              SSData[:,SSFields.index(parm)]) / float(scale)
#                    ax.plot(U_Refs,y_SS,'k:')
#                
#                # plot data
#                for iI in range(len(WindParms[1])):
#                    I = WindParms[1][iI]
#                    x = np.empty((len(UPlot),4))
#                    x[:,0],x[:,1],x[:,2],x[:,3] = UPlot,I,logL,rho
#                    Xv = jr.myvander(x,ps)
#                    yplot = np.dot(Xv,cs)
#                    
#                    ax.plot(UPlot,yplot,
#                            label='$I$ = {:.1f}'.format(WindParms[1][iI]))
#                    
#                # prettify axes
#                ax.set_xlim(ULim)
#                plt.locator_params(axis='both',nbins=TickBins,tight=True)
#                jr.removeSpines(ax)
#                
#                # put x and y labels on left column and bottom row only
#                if iRho < len(WindParms[3])-1:
#                    ax.set_xticklabels([])
#                if iL > 0:
#                    ax.set_yticklabels([])
#                    
#                
#                # create legend
#                if (iL == 0) and (iRho == 0):
#                    plt.legend(bbox_to_anchor=LegLoc,
#                               bbox_transform=plt.gcf().transFigure, loc=3,
#                               ncol=5, mode='expand', borderaxespad=0.,
#                               fontsize=LegFS)
#        
#        # scale subplots and add text labels
#        xbord,ybord = 0.07,0.025
#        plt.tight_layout(rect=[xbord,ybord,1.01,0.94])
#        plt.figtext(0.5+xbord/2.,0.02,'Mean Wind Speed [m/s]',
#                    ha='center',va='center',fontsize='small')
#        plt.figtext(0.08,0.5+ybord/2.,'{:s} {:s} [{:s}]'.format(stat,parm,units),
#                    ha='center',va='center',rotation='vertical',fontsize='small')
#        
#        for iRho in range(len(WindParms[3])):
#            plt.figtext(0.025,0.85-0.88*iRho/len(WindParms[3]),
#                        r'$\rho$ = {:.1f}'.format(WindParms[3][iRho]),
#                        ha='center',va='center',rotation='vertical',fontsize='small')
#        for iL in range(len(WindParms[2])):
#            plt.figtext(0.23 + 0.88*iL/len(WindParms[2]),0.98,
#                        r'log$_{10}$L = ' + '{:.1f}'.format(WindParms[2][iL]),
#                        ha='center',va='center',fontsize='small')
#                    
#        if savefig:
#            FigTitle = 'plot-RSM_lines_vsI_{:s}_{:s}_{:s}.eps'.format(TurbName,parm,stat)
#            SavePath = os.path.join(SaveDir,FigTitle)
#            fig.savefig(SavePath)
#            print('Figure {:s} saved.'.format(FigTitle))
#            
#        # ================= plot data vs rho =====================
#        
#        # initialize figure
#        fig = plt.figure(FigNum+1,figsize=FigSize)
#        plt.clf()
#        
#        # plot data
#        for iI in range(len(WindParms[1])):
#            I = WindParms[1][iI]
#            
#            for iL in range(len(WindParms[2])):
#                logL = WindParms[2][iL]
#                    
#                # create axes
#                iplot = iI*len(WindParms[2]) + iL + 1
#                ax = fig.add_subplot(len(WindParms[1]),len(WindParms[2]),iplot)
#                
#                # plot steady-state if applicable
#                if stat in ['mean','max']:
#                    U_Refs = np.linspace(min(WindParms[0]),max(WindParms[0]),
#                                         20)
#                    U_hubs = U_Refs*(HubHeight/zRef)**shear
#                    y_SS   = np.interp(U_hubs,SSData[:,SSFields.index('WindVxi')],
#                                              SSData[:,SSFields.index(parm)]) / float(scale)
#                    ax.plot(U_Refs,y_SS,'k:')
#                
#                # plot data
#                for iRho in range(len(WindParms[3])):
#                    rho = WindParms[3][iRho]
#                    x = np.empty((len(UPlot),4))
#                    x[:,0],x[:,1],x[:,2],x[:,3] = UPlot,I,logL,rho
#                    Xv = jr.myvander(x,ps)
#                    yplot = np.dot(Xv,cs)
#                    
#                    ax.plot(UPlot,yplot,
#                            label=r'$\rho$ = {:.1f}'.format(WindParms[3][iRho]))
#                    
#                # prettify axes
#                ax.set_xlim(ULim)
#                plt.locator_params(axis='both',nbins=TickBins,tight=True)
#                jr.removeSpines(ax)
#                
#                # put x and y labels on left column and bottom row only
#                if iI < len(WindParms[1])-1:
#                    ax.set_xticklabels([])
#                if iL > 0:
#                    ax.set_yticklabels([])
#                
#                # create legend
#                if (iL == 0) and (iI == 0):
#                    plt.legend(bbox_to_anchor=LegLoc,
#                               bbox_transform=plt.gcf().transFigure, loc=3,
#                               ncol=5, mode='expand', borderaxespad=0.,
#                               fontsize=LegFS)
#        
#        # scale subplots and add text labels
#        xbord,ybord = 0.07,0.025
#        plt.tight_layout(rect=[xbord,ybord,1.01,0.94])
#        plt.figtext(0.5+xbord/2.,0.02,'Mean Wind Speed [m/s]',
#                    ha='center',va='center',fontsize='small')
#        plt.figtext(0.08,0.5+ybord/2.,'{:s} {:s} [{:s}]'.format(stat,parm,units),
#                    ha='center',va='center',rotation='vertical',fontsize='small')
#        
#        for iI in range(len(WindParms[1])):
#            plt.figtext(0.025,0.85-0.88*iI/len(WindParms[1]),
#                        r'$I$ = {:.1f}'.format(WindParms[1][iI]),
#                        ha='center',va='center',rotation='vertical',fontsize='small')
#        for iL in range(len(WindParms[2])):
#            plt.figtext(0.23 + 0.88*iL/len(WindParms[2]),0.98,
#                        r'log$_{10}$L = ' + '{:.1f}'.format(WindParms[2][iL]),
#                        ha='center',va='center',fontsize='small')
#                           
#        if savefig:
#            FigTitle = 'plot-RSM_lines_vsRho_{:s}_{:s}_{:s}.eps'.format(TurbName,parm,stat)
#            SavePath = os.path.join(SaveDir,FigTitle)
#            fig.savefig(SavePath)
#            print('Figure {:s} saved.'.format(FigTitle))
#            
#            
#        # ================= plot data vs L =====================
#        
#        # initialize figure
#        fig = plt.figure(FigNum+2,figsize=FigSize)
#        plt.clf()
#        
#        # plot data
#        for iI in range(len(WindParms[1])):
#            I = WindParms[1][iI]
#            
#            for iRho in range(len(WindParms[3])):
#                rho = WindParms[3][iRho]
#        
#                # create axes
#                iplot = iI*len(WindParms[3]) + iRho + 1
#                ax = fig.add_subplot(len(WindParms[1]),len(WindParms[3]),iplot)
#                
#                # plot steady-state if applicable
#                if stat in ['mean','max']:
#                    U_Refs = np.linspace(min(WindParms[0]),max(WindParms[0]),
#                                         20)
#                    U_hubs = U_Refs*(HubHeight/zRef)**shear
#                    y_SS   = np.interp(U_hubs,SSData[:,SSFields.index('WindVxi')],
#                                              SSData[:,SSFields.index(parm)]) / float(scale)
#                    ax.plot(U_Refs,y_SS,'k:')
#                
#                # plot data
#                for iL in range(len(WindParms[2])):
#                    logL = WindParms[2][iL]
#                    x = np.empty((len(UPlot),4))
#                    x[:,0],x[:,1],x[:,2],x[:,3] = UPlot,I,logL,rho
#                    Xv = jr.myvander(x,ps)
#                    yplot = np.dot(Xv,cs)
#                    
#                    ax.plot(UPlot,yplot,
#                            label='$L$ = $10^{' + '{:.1f}'.format(WindParms[2][iL]) + \
#                                '}$')
#                    
#                # prettify axes
#                ax.set_xlim(ULim)
#                plt.locator_params(axis='both',nbins=TickBins,tight=True)
#                jr.removeSpines(ax)
#                
#                # put x and y labels on left column and bottom row only
#                if iI < len(WindParms[1])-1:
#                    ax.set_xticklabels([])
#                if iRho > 0:
#                    ax.set_yticklabels([])
#                
#                # create legend
#                if (iRho == 0) and (iI == 0):
#                    plt.legend(bbox_to_anchor=LegLoc,
#                               bbox_transform=plt.gcf().transFigure, loc=3,
#                               ncol=5, mode='expand', borderaxespad=0.,
#                               fontsize=LegFS)
#        
#        # scale subplots and add text labels
#        xbord,ybord = 0.07,0.025
#        plt.tight_layout(rect=[xbord,ybord,1.01,0.94])
#        plt.figtext(0.5+xbord/2.,0.02,'Mean Wind Speed [m/s]',
#                    ha='center',va='center',fontsize='small')
#        plt.figtext(0.08,0.5+ybord/2.,'{:s} {:s} [{:s}]'.format(stat,parm,units),
#                    ha='center',va='center',rotation='vertical',fontsize='small')
#        
#        for iI in range(len(WindParms[1])):
#            plt.figtext(0.025,0.85-0.88*iI/len(WindParms[1]),
#                        r'$I$ = {:.1f}'.format(WindParms[1][iI]),
#                        ha='center',va='center',rotation='vertical',fontsize='small')
#        for iRho in range(len(WindParms[2])):
#            plt.figtext(0.23 + 0.88*iRho/len(WindParms[2]),0.98,
#                        '$\\rho$ = ' + '{:.1f}'.format(WindParms[3][iRho]),
#                        ha='center',va='center',fontsize='small')
#        
#        if savefig:
#            FigTitle = 'plot-RSM_lines_vsL_{:s}_{:s}_{:s}.eps'.format(TurbName,parm,stat)
#            SavePath = os.path.join(SaveDir,FigTitle)
#            fig.savefig(SavePath)
#            print('Figure {:s} saved.'.format(FigTitle))
            
