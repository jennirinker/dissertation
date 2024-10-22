"""
Plot statistic versus rho for different L and sigma
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
import os
import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm


# define turbine name and run name
#TurbNames = ['WP0.75A08V00','WP1.5A08V03',
#              'WP3.0A02V02','WP5.0A04V00']
#RunNames = ['Peregrine','TestRun','SmallRun','BigRun2']
TurbName = 'WP5.0A04V00'
RunName  = 'BigRun2'

FigNum  = 6
FigSize = (6.0,7.5)
alpha   = 0.05

# statistic and value to fit RSM to
#stat,parm,units,scale = 'max','RootMFlp1','MN-m',1000.
stat,parm,units,scale = 'DEL-h','RootMFlp1','MN-m',1000.
#stat,parm,units,scale = 'max','HSShftTq','kN-m',1
#stat,parm,units,scale = 'DEL-h','HSShftTq','kN-m',1
#stat,parm,units,scale = 'max','TwrBsMyt','MN-m',1000
#stat,parm,units,scale = 'DEL-h','TwrBsMyt','MN-m',1000

# base directory where the stats are stored
BaseStatDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                'processed_data\\proc_stats'


# -----------------------------------------------------------------------------

def StatsErr(x,y,p_i,
             alpha=0.05):
    """ Sum-of-squared-error for data x and y with polynomial orders p_i
    """
    
    # perform OLS for all data
    ps_all = jr.GetAllPowers(p_i)
    Xv_all = jr.myvander(x,ps_all)
    results = sm.OLS(y, Xv_all).fit()
    
    # extract significant coefficients
    ps_red = ps_all[results.pvalues <= alpha]
    Xv_red = jr.myvander(x,ps_red)
    
    # fit OLS with reduced coefficients
    cs_red = sm.OLS(y, Xv_red).fit().params
        
    # get reduced model
    yhat   = np.dot(Xv_red,cs_red)
    
    # get residuals and sum
    e    = y - yhat
    err  = np.mean(e ** 2)
    
    return err

# get wind parameters for that run
WindParms = jr.RunName2WindParms(RunName)
URefs, Is, Ls, rhos, n_dups = WindParms['URefs'],WindParms['Is'], \
                              WindParms['Ls'],WindParms['rhos'], \
                              WindParms['n_dups']
WindParms = [URefs,Is,np.log10(Ls),rhos]
WindParmStr = ['$U$','$\sigma_u$','log$_{10}(L)$',r'$\rho$']

# load the stats data
x, y = jr.LoadFASTStats(RunName,TurbName,stat,parm)

# ================= fit polynomial surface =====================

# parameterize function for data
ErrFunc = lambda p: StatsErr(x,y,p,alpha=alpha)

p0 = np.zeros(x.shape[1])
results = jr.DiscreteOpt(ErrFunc,p0,
                         verbose=1)
p_i = results['p_out']
print(p_i)
#p_i = [6,2,2,2]
#p_i = [10,3,5,2]

# get significant powers and coefficients
ps_all = jr.GetAllPowers(p_i)
Xv_all = jr.myvander(x,ps_all)
results = sm.OLS(y, Xv_all).fit()
cs_all  = results.params
ps_red = ps_all[results.pvalues <= alpha]
Xv_red = jr.myvander(x,ps_red)
cs_red = sm.OLS(y, Xv_red).fit().params

yhat   = np.dot(Xv_red,cs_red)
es_red = y - yhat



# **********************************************************
# scale data
e = es_red/scale

# ================= plot data vs I =====================

# initialize figure
fig = plt.figure(FigNum,figsize=FigSize)
plt.clf()
        
# plot data
emax = max(-np.floor(e.min()), np.ceil(e.max()))
for iL in range(len(WindParms[2])):
    logL = WindParms[2][iL]
    for iRho in range(len(WindParms[3])):
        rho = WindParms[3][iRho]
        
        X,Y  = np.meshgrid(WindParms[0],WindParms[1])
        Z    = np.empty(X.shape)
        Zerr = np.empty(X.shape)
        
        for iU, iI in [(a,b) for a in range(len(WindParms[0])) \
                             for b in range(len(WindParms[1]))]:
            
            U, I = WindParms[0][iU], WindParms[1][iI]
            mask = np.logical_and(np.logical_and(np.logical_and(x[:,0]==U,
                                                        x[:,1]==I),
                                         x[:,2]==logL),
                          x[:,3]==rho)
            e_data = e[mask]
            Z[iI,iU]    = np.mean(e_data)
            Zerr[iI,iU] = np.std(e_data)
            
        # create axes
        iplot = iRho*len(WindParms[2]) + iL + 1
        ax = fig.add_subplot(len(WindParms[3]),len(WindParms[2]),iplot)
        
        # plot data
        for iI in range(len(WindParms[1])):
            ax.errorbar(WindParms[0],Z[iI,:],yerr=Zerr[iI,:],
                        label='$I$ = {:.1f}'.format(WindParms[1][iI]))
            
        # prettify axes
        ax.set_ylim([-emax,emax])
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
               
    
# ================= plot data vs rho =====================

# initialize figure
fig = plt.figure(FigNum+1,figsize=FigSize)
plt.clf()

# plot data
#vmin, vmax = np.floor(y.min()/100.)*100., np.ceil(y.max()/100.)*100
emax = max(-np.floor(e.min()), np.ceil(e.max()))
for iL in range(len(WindParms[2])):
    logL = WindParms[2][iL]
    for iI in range(len(WindParms[1])):
        I = WindParms[1][iI]
        
        X,Y  = np.meshgrid(WindParms[0],WindParms[3])
        Z    = np.empty(X.shape)
        Zerr = np.empty(X.shape)
        
        for iU, iRho in [(a,b) for a in range(len(WindParms[0])) \
                             for b in range(len(WindParms[3]))]:
            
            U, rho = WindParms[0][iU], WindParms[3][iRho]
            mask = np.logical_and(np.logical_and(np.logical_and(x[:,0]==U,
                                                        x[:,1]==I),
                                         x[:,2]==logL),
                          x[:,3]==rho)
            e_data = e[mask]
            Z[iRho,iU]    = np.mean(e_data)
            Zerr[iRho,iU] = np.std(e_data)
            
        # create axes
        iplot = iI*len(WindParms[2]) + iL + 1
        ax = fig.add_subplot(len(WindParms[1]),len(WindParms[2]),iplot)
        
        # plot data
        for iRho in range(len(WindParms[3])):
            ax.errorbar(WindParms[0],Z[iRho,:],yerr=Zerr[iRho,:],
                        label=r'$\rho$ = {:.1f}'.format(WindParms[3][iRho]))
            
        # prettify axes
        ax.set_ylim([-emax,emax])
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





