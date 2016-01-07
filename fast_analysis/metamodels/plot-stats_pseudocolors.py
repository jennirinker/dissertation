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


# define turbine name and run name
#TurbNames = ['WP0.75A08V00','WP1.5A08V03',
#              'WP3.0A02V02','WP5.0A04V00']
#RunNames = ['Peregrine','TestRun','SmallRun','BigRun2']
TurbName = 'WP5.0A04V00'
RunName  = 'BigRun2'

FigNum = 1

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
stat,parm = 'max','RootMFlp1'
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

# ================= plot data =====================

## get mean of each set of parameters
#n_pts = np.prod(np.array([len(l) for l in WindParms]))
#x_mean = np.empty((n_pts,x.shape[1]))
#y_mean = np.empty(n_pts)
#y_std  = np.empty(n_pts)
#i_ct   = 0
#for U,I,logL,rho in [(a,b,c,d) for a in WindParms[0] for b in WindParms[1] \
#                               for c in WindParms[2] for d in WindParms[3]]:
#    mask = np.logical_and(np.logical_and(np.logical_and(x[:,0]==U,
#                                                        x[:,1]==I),
#                                         x[:,2]==logL),
#                          x[:,3]==rho)
#    x_mean[i_ct] = [U,I,logL,rho]
#    y_mean[i_ct] = np.mean(y[mask])
#    y_std[i_ct]  = np.std(y[mask])
#    
#    i_ct += 1
                

# initialize figure
#fig = plt.figure(FigNum,figsize=(6,6))
fig = plt.figure(FigNum,figsize=(10,10))
plt.clf()

# plot data
vmin, vmax = np.floor(y.min()/100.)*100., np.ceil(y.max()/100.)*100
levels = np.linspace(vmin,vmax,10)
for iL in range(len(WindParms[2])):
    logL = WindParms[2][iL]
    for iRho in range(len(WindParms[3])):
        rho = WindParms[3][iRho]
        
        X,Y = np.meshgrid(WindParms[0],WindParms[1])
        Z   = np.empty(X.shape)
        
        for iU, iI in [(a,b) for a in range(len(WindParms[0])) \
                             for b in range(len(WindParms[1]))]:
            
            U, I = WindParms[0][iU], WindParms[1][iI]
            mask = np.logical_and(np.logical_and(np.logical_and(x[:,0]==U,
                                                        x[:,1]==I),
                                         x[:,2]==logL),
                          x[:,3]==rho)
            y_data = y[mask]
            Z[iI,iU] = np.mean(y_data)
            
        iplot = iRho*len(WindParms[2]) + iL + 1
        ax = fig.add_subplot(len(WindParms[3]),len(WindParms[2]),iplot)
        cnt = ax.contourf(X,Y,Z,cmap='Reds',levels=levels)
        cbar = plt.colorbar(cnt)
#        cbar.set_clim([vmin,vmax])

plt.tight_layout()