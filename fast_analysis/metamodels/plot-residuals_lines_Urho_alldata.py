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
import matplotlib.ticker as mtick

# Duke style
plt.style.use(jr.stylepath('duke_paper'))

# define turbine name and run name
#TurbNames = ['WP0.75A08V00','WP1.5A08V03',
#              'WP3.0A02V02','WP5.0A04V00']
#RunNames = ['Peregrine','TestRun','SmallRun','BigRun2']
TurbName = 'WP5.0A04V00'
RunName  = 'BigRun2'

FigNum = 4

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
WindParmStr = ['$U$','$\sigma_u$','log$_{10}(L)$',r'$\rho$']

# load the stats data
DictPath   = os.path.join(BaseStatDir,RunName,TurbName + '_stats.mat')
stats_dict = scio.loadmat(DictPath,squeeze_me=True)
proc_stats = stats_dict['proc_stats']
calc_stats = [s.rstrip() for s in stats_dict['calc_stats']]
fnames     = [s.rstrip() for s in stats_dict['fnames']]
fields     = [s.rstrip() for s in stats_dict['fields']]
n_fields = len(fields)

# statistic and value to fit RSM to
stat,parm,units,scale = 'max','RootMFlp1','MN-m',1000.
#stat,parm = 'max','RootMEdg1'
#stat,parm = 'max','RotTorq'
#stat,parm = 'max','HSShftTq'
#stat,parm = 'max','YawBrMzp'
#stat,parm = 'max','TwrBsMyt'

# extract data
y = proc_stats[:,calc_stats.index(stat)*n_fields + \
                 fields.index(parm)]                # output data
x = np.empty((y.size,4))
for i_f in range(y.size):
    file_id =  fnames[i_f].rstrip('.out').split('_')[1]
    for i_p in range(len(WindParms)):
        x[i_f,i_p] = WindParms[i_p][int(file_id[i_p],16)]   # hex to int


# calculate fit
p_i = [7,2,2,2]
betas, ps = jr.polyregression(x,y,p_i)
X    = jr.myvander(x,ps)
yhat = np.dot(X,betas)
e = y - yhat

# **********************************************************
# scale data
e = e/scale

# **********************************************************

# ================= plot data =====================

# initialize figure
fig = plt.figure(FigNum,figsize=(6,6))
plt.clf()

# initialize figure
#fig = plt.figure(FigNum,figsize=(6,6))
fig = plt.figure(FigNum,figsize=(10,10))
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
    

