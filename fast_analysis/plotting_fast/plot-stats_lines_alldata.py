"""
Plot statistic versus rho for different L and sigma
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
import os, json
import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt

# Duke style
plt.style.use(jr.stylepath('duke_paper'))

# define turbine name and run name
#TurbNames = ['WP0.75A08V00','WP1.5A08V03',
#              'WP3.0A02V02','WP5.0A04V00']
#RunNames = ['Peregrine','TestRun','SmallRun','BigRun2']
TurbName = 'WP3.0A02V02'
RunName  = 'BigRun2'

FigNum = 3
#FigSize = (6.0,8.5)    # fulll page
FigSize = (4.25,6.0)    # half page
LegLoc = (0.03, 0.92, 0.96, .05)
LegFS = 'x-small'
ULim = [4,24]
TickBins = 4

zRef  = 90.
shear = 0.2

# base directory where the stats are stored
BaseStatDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                'processed_data\\proc_stats'
BaseTurbDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
                'dissertation\\FAST_models\\FAST7'


# -----------------------------------------------------------------------------

# get wind parameters for that run
WindParms = jr.RunName2WindParms(RunName)
URefs, Is, Ls, rhos, n_dups = WindParms['URefs'],WindParms['Is'], \
                              WindParms['Ls'],WindParms['rhos'], \
                              WindParms['n_dups']
WindParms = [URefs,Is,np.log10(Ls),rhos]
WindParmStr = ['$U$','$\sigma_u$','log$_{10}(L)$',r'$\rho$']

# statistic and value to fit RSM to
#stat,parm,units,scale = 'max','RootMFlp1','MN-m',1000.
#stat,parm = 'max','RootMEdg1'
#stat,parm = 'max','RotTorq'
#stat,parm = 'max','HSShftTq'
#stat,parm = 'max','YawBrMzp'
#stat,parm = 'max','TwrBsMyt'
stat,parm,units,scale = 'DEL-m','RootMFlp1','MN-m',1000.

x, y = jr.LoadFASTStats(RunName,TurbName,stat,parm,
                        scale=scale)

# get unique means and standard deviations
x_uniq, err_mean, err_std = jr.FASTUniqueStats(x,y,RunName)

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

# ================= plot data vs I =====================

# initialize figure
fig = plt.figure(FigNum,figsize=FigSize)
plt.clf()

# plot data
for iRho in range(len(WindParms[3])):
    rho = WindParms[3][iRho]
    
    # mask row data to determine axes limits
    emean_row = err_mean[x_uniq[:,3] == rho]
    estd_row  = err_std[x_uniq[:,3] == rho]
    emax = emean_row + estd_row
    vmin, vmax = np.floor(emax.min()), np.ceil(emax.max())
    
    for iL in range(len(WindParms[2])):
        logL = WindParms[2][iL]
        
        X,Y  = np.meshgrid(WindParms[0],WindParms[1])
        Z    = np.empty(X.shape)
        Zerr = np.empty(X.shape)
        
        for iU, iI in [(a,b) for a in range(len(WindParms[0])) \
                             for b in range(len(WindParms[1]))]:
            
            U, I = WindParms[0][iU], WindParms[1][iI]
            mask = np.logical_and(np.logical_and(np.logical_and(x_uniq[:,0]==U,
                                                        x_uniq[:,1]==I),
                                         x_uniq[:,2]==logL),
                          x_uniq[:,3]==rho)
            Z[iI,iU]    = err_mean[mask]
            Zerr[iI,iU] = err_std[mask]
            
        # create axes
        iplot = iRho*len(WindParms[2]) + iL + 1
        ax = fig.add_subplot(len(WindParms[3]),len(WindParms[2]),iplot)
                
        # plot steady-state if applicable
        if stat in ['mean','max']:
            U_Refs = np.linspace(min(WindParms[0]),max(WindParms[0]),
                                 20)
            U_hubs = U_Refs*(HubHeight/zRef)**shear
            y_SS   = np.interp(U_hubs,SSData[:,SSFields.index('WindVxi')],
                                      SSData[:,SSFields.index(parm)]) / float(scale)
            ax.plot(U_Refs,y_SS,'k:')
        
        # plot data
        for iI in range(len(WindParms[1])):
            ax.errorbar(WindParms[0],Z[iI,:],yerr=Zerr[iI,:],
                        label='$I$ = {:.1f}'.format(WindParms[1][iI]))
            
        # prettify axes
        ax.set_ylim([vmin,vmax])
        ax.set_xlim(ULim)
        plt.locator_params(axis='both',nbins=TickBins,tight=True)
        jr.removeSpines(ax)
        
        # put x and y labels on left column and bottom row only
        if iRho < len(WindParms[3])-1:
            ax.set_xticklabels([])
        if iL > 0:
            ax.set_yticklabels([])
            
        
        # create legend
        if (iL == 0) and (iRho == 0):
            plt.legend(bbox_to_anchor=LegLoc,
                       bbox_transform=plt.gcf().transFigure, loc=3,
                       ncol=5, mode='expand', borderaxespad=0.,
                       fontsize=LegFS)

# scale subplots and add text labels
xbord,ybord = 0.07,0.025
plt.tight_layout(rect=[xbord,ybord,1.01,0.94])
plt.figtext(0.5+xbord/2.,0.02,'Mean Wind Speed [m/s]',
            ha='center',va='center',fontsize='small')
plt.figtext(0.07,0.5+ybord/2.,'{:s} {:s} [{:s}]'.format(stat,parm,units),
            ha='center',va='center',rotation='vertical',fontsize='small')

for iRho in range(len(WindParms[3])):
    plt.figtext(0.025,0.85-0.88*iRho/len(WindParms[3]),
                r'$\rho$ = {:.1f}'.format(WindParms[3][iRho]),
                ha='center',va='center',rotation='vertical',fontsize='small')
for iL in range(len(WindParms[2])):
    plt.figtext(0.23 + 0.88*iL/len(WindParms[2]),0.98,
                r'log$_{10}$L = ' + '{:.1f}'.format(WindParms[2][iL]),
                ha='center',va='center',fontsize='small')
    

# ================= plot data vs rho =====================

# initialize figure
fig = plt.figure(FigNum+1,figsize=FigSize)
plt.clf()

# plot data
for iI in range(len(WindParms[1])):
    I = WindParms[1][iI]
    
    # mask row data to determine axes limits
    emean_row = err_mean[x_uniq[:,1] == I]
    estd_row  = err_std[x_uniq[:,1] == I]
    emax = emean_row + estd_row
    vmin, vmax = np.floor(emax.min()), np.ceil(emax.max())
    
    for iL in range(len(WindParms[2])):
        logL = WindParms[2][iL]
        
        X,Y  = np.meshgrid(WindParms[0],WindParms[3])
        Z    = np.empty(X.shape)
        Zerr = np.empty(X.shape)
        
        for iU, iRho in [(a,b) for a in range(len(WindParms[0])) \
                             for b in range(len(WindParms[3]))]:
            
            U, rho = WindParms[0][iU], WindParms[3][iRho]
            mask = np.logical_and(np.logical_and(np.logical_and(x_uniq[:,0]==U,
                                                        x_uniq[:,1]==I),
                                         x_uniq[:,2]==logL),
                          x_uniq[:,3]==rho)
            Z[iRho,iU]    = err_mean[mask]
            Zerr[iRho,iU] = err_std[mask]
            
        # create axes
        iplot = iI*len(WindParms[2]) + iL + 1
        ax = fig.add_subplot(len(WindParms[1]),len(WindParms[2]),iplot)
                
        # plot steady-state if applicable
        if stat in ['mean','max']:
            U_Refs = np.linspace(min(WindParms[0]),max(WindParms[0]),
                                 20)
            U_hubs = U_Refs*(HubHeight/zRef)**shear
            y_SS   = np.interp(U_hubs,SSData[:,SSFields.index('WindVxi')],
                                      SSData[:,SSFields.index(parm)]) / float(scale)
            ax.plot(U_Refs,y_SS,'k:')
        
        # plot data
        for iRho in range(len(WindParms[3])):
            ax.errorbar(WindParms[0],Z[iRho,:],yerr=Zerr[iRho,:],
                        label=r'$\rho$ = {:.1f}'.format(WindParms[3][iRho]))
            
        # prettify axes
        ax.set_ylim([vmin,vmax])
        ax.set_xlim(ULim)
        plt.locator_params(axis='both',nbins=TickBins,tight=True)
        jr.removeSpines(ax)
        
        # put x and y labels on left column and bottom row only
        if iI < len(WindParms[1])-1:
            ax.set_xticklabels([])
        if iL > 0:
            ax.set_yticklabels([])
        
        # create legend
        if (iL == 0) and (iI == 0):
            plt.legend(bbox_to_anchor=LegLoc,
                       bbox_transform=plt.gcf().transFigure, loc=3,
                       ncol=5, mode='expand', borderaxespad=0.,
                       fontsize=LegFS)

# scale subplots and add text labels
xbord,ybord = 0.07,0.025
plt.tight_layout(rect=[xbord,ybord,1.01,0.94])
plt.figtext(0.5+xbord/2.,0.02,'Mean Wind Speed [m/s]',
            ha='center',va='center',fontsize='small')
plt.figtext(0.07,0.5+ybord/2.,'{:s} {:s} [{:s}]'.format(stat,parm,units),
            ha='center',va='center',rotation='vertical',fontsize='small')

for iI in range(len(WindParms[1])):
    plt.figtext(0.025,0.85-0.88*iI/len(WindParms[1]),
                r'$I$ = {:.1f}'.format(WindParms[1][iI]),
                ha='center',va='center',rotation='vertical',fontsize='small')
for iL in range(len(WindParms[2])):
    plt.figtext(0.23 + 0.88*iL/len(WindParms[2]),0.98,
                r'log$_{10}$L = ' + '{:.1f}'.format(WindParms[2][iL]),
                ha='center',va='center',fontsize='small')
           
# ================= plot data vs L =====================

# initialize figure
fig = plt.figure(FigNum+2,figsize=FigSize)
plt.clf()

# plot data
for iI in range(len(WindParms[1])):
    I = WindParms[1][iI]
    
    # mask row data to determine axes limits
    emean_row = err_mean[x_uniq[:,1] == I]
    estd_row  = err_std[x_uniq[:,1] == I]
    emax = emean_row + estd_row
    vmin, vmax = np.floor(emax.min()), np.ceil(emax.max())
    
    for iRho in range(len(WindParms[3])):
        rho = WindParms[3][iRho]

        
        X,Y  = np.meshgrid(WindParms[0],WindParms[2])
        Z    = np.empty(X.shape)
        Zerr = np.empty(X.shape)
        
        for iU, iL in [(a,b) for a in range(len(WindParms[0])) \
                             for b in range(len(WindParms[2]))]:
            
            U, logL = WindParms[0][iU], WindParms[2][iL]
            mask = np.logical_and(np.logical_and(np.logical_and(x_uniq[:,0]==U,
                                                        x_uniq[:,1]==I),
                                         x_uniq[:,2]==logL),
                          x_uniq[:,3]==rho)
            Z[iL,iU]    = err_mean[mask]
            Zerr[iL,iU] = err_std[mask]
            
        # create axes
        iplot = iI*len(WindParms[3]) + iRho + 1
        ax = fig.add_subplot(len(WindParms[1]),len(WindParms[3]),iplot)
                
        # plot steady-state if applicable
        if stat in ['mean','max']:
            U_Refs = np.linspace(min(WindParms[0]),max(WindParms[0]),
                                 20)
            U_hubs = U_Refs*(HubHeight/zRef)**shear
            y_SS   = np.interp(U_hubs,SSData[:,SSFields.index('WindVxi')],
                                      SSData[:,SSFields.index(parm)]) / float(scale)
            ax.plot(U_Refs,y_SS,'k:')
        
        # plot data
        for iL in range(len(WindParms[2])):
            ax.errorbar(WindParms[0],Z[iL,:],yerr=Zerr[iL,:],
                        label='$L$ = $10^{' + '{:.1f}'.format(WindParms[2][iL]) + \
                            '}$')
            
        # prettify axes
        ax.set_ylim([vmin,vmax])
        ax.set_xlim(ULim)
        plt.locator_params(axis='both',nbins=TickBins,tight=True)
        jr.removeSpines(ax)
        
        # put x and y labels on left column and bottom row only
        if iI < len(WindParms[1])-1:
            ax.set_xticklabels([])
        if iRho > 0:
            ax.set_yticklabels([])
        
        # create legend
        if (iRho == 0) and (iI == 0):
            plt.legend(bbox_to_anchor=LegLoc,
                       bbox_transform=plt.gcf().transFigure, loc=3,
                       ncol=5, mode='expand', borderaxespad=0.,
                       fontsize=LegFS)

# scale subplots and add text labels
xbord,ybord = 0.07,0.025
plt.tight_layout(rect=[xbord,ybord,1.01,0.94])
plt.figtext(0.5+xbord/2.,0.02,'Mean Wind Speed [m/s]',
            ha='center',va='center',fontsize='small')
plt.figtext(0.07,0.5+ybord/2.,'{:s} {:s} [{:s}]'.format(stat,parm,units),
            ha='center',va='center',rotation='vertical',fontsize='small')

for iI in range(len(WindParms[1])):
    plt.figtext(0.025,0.85-0.88*iI/len(WindParms[1]),
                r'$I$ = {:.1f}'.format(WindParms[1][iI]),
                ha='center',va='center',rotation='vertical',fontsize='small')
for iRho in range(len(WindParms[2])):
    plt.figtext(0.23 + 0.88*iRho/len(WindParms[2]),0.98,
                '$\\rho$ = ' + '{:.1f}'.format(WindParms[3][iRho]),
                ha='center',va='center',fontsize='small')
