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

# initialize figure
fig = plt.figure(FigNum,figsize=(6,6))
plt.clf()

# plotting idxs
ips1,ips2 = 0, 1            # indices of screening parameters
ip   = 3                    # indices of plotting parameter
il   = 2                    # indices of line parameter
is1s = [0,-1]               # indices for SP1
is2s = [0,-1]               # indices for SP2

# loop through U and sigma
n1, n2 = len(is1s), len(is2s)
xpstr = WindParmStr[ip]
for i in range(n1):
    is1   = is1s[i]
    x1    = WindParms[ips1][is1]
    x1str = WindParmStr[ips1]
    for j in range(n2):
        is2   = is2s[j]
        x2    = WindParms[ips2][is2]
        x2str = WindParmStr[ips2]
        
        mask   = np.logical_and(x[:,ips1]==x1,x[:,ips2]==x2)
        x_mask = x[mask,:]
        y_mask = y[mask]
        
        # create axes
        i_sp = j*n2+i
        ax   = fig.add_subplot(n2,n1,i_sp+1)
        
        # loop through line parameters
        lps = WindParms[il]
        for i_lp in range(len(lps)):
            lp = lps[i_lp]
            
            # get line data
            idcs_data = x_mask[:,il] == lp
            x_data    = x_mask[idcs_data,ip]
            y_data    = y_mask[idcs_data]
            
            # take average of multiple values
            x_plot = np.unique(x_data)
            y_plot = np.empty(x_plot.shape)
            y_errb = np.empty(x_plot.shape)
            for i_xuniq in range(len(x_plot)):
                x_uniq = x_plot[i_xuniq]
                y_plot[i_xuniq] = np.mean(y_data[x_data == x_uniq])
                y_errb[i_xuniq] = np.std(y_data[x_data == x_uniq])
            
#            ax.errorbar(x_plot,y_plot,yerr=y_errb)
            ax.errorbar(x_plot+0.005*i_lp,y_plot,yerr=y_errb,
                        label='{:.1f}'.format(lp))
                        
        # prettify axes
        ax.set_title('{:s} = {:.1f}, {:s} = {:.1f}'.format(x1str,x1,x2str,x2))
        ax.set_xlabel('{:s}'.format(xpstr))
        ax.set_xlim([-0.02,0.45])
#        ax.set_ylim(ylim)
        
        if i_sp == n1-1:
            plt.legend(fontsize='x-small',loc=0)

plt.tight_layout()


# ================= plot data =====================

# initialize figure
fig = plt.figure(FigNum+1,figsize=(6,6))
plt.clf()

# plotting idxs
ips1,ips2 = 2, 3            # indices of screening parameters
ip   = 0                    # indices of plotting parameter
il   = 1                    # indices of line parameter
is1s = [0,-1]               # indices for SP1
is2s = [0,-1]               # indices for SP2

# loop through U and sigma
n1, n2 = len(is1s), len(is2s)
xpstr = WindParmStr[ip]
for i in range(n1):
    is1   = is1s[i]
    x1    = WindParms[ips1][is1]
    x1str = WindParmStr[ips1]
    for j in range(n2):
        is2   = is2s[j]
        x2    = WindParms[ips2][is2]
        x2str = WindParmStr[ips2]
        
        mask   = np.logical_and(x[:,ips1]==x1,x[:,ips2]==x2)
        x_mask = x[mask,:]
        y_mask = y[mask]
        
        # create axes
        i_sp = j*n2+i
        ax   = fig.add_subplot(n2,n1,i_sp+1)
        
        # loop through line parameters
        lps = WindParms[il]
        for i_lp in range(len(lps)):
            lp = lps[i_lp]
            
            # get line data
            idcs_data = x_mask[:,il] == lp
            x_data    = x_mask[idcs_data,ip]
            y_data    = y_mask[idcs_data]
            
            # take average of multiple values
            x_plot = np.unique(x_data)
            y_plot = np.empty(x_plot.shape)
            y_errb = np.empty(x_plot.shape)
            for i_xuniq in range(len(x_plot)):
                x_uniq = x_plot[i_xuniq]
                y_plot[i_xuniq] = np.mean(y_data[x_data == x_uniq])
                y_errb[i_xuniq] = np.std(y_data[x_data == x_uniq])
            
#            ax.errorbar(x_plot,y_plot,yerr=y_errb)
            ax.errorbar(x_plot,y_plot,yerr=y_errb,
                        label='{:.1f}'.format(lp))
        
        # prettify axes
        ax.set_title('{:s} = {:.1f}, {:s} = {:.1f}'.format(x1str,x1,x2str,x2))
        ax.set_xlabel('{:s}'.format(xpstr))
#        ax.set_ylim(ylim)
        
        if i_sp == n1-1:
            plt.legend(fontsize='x-small',loc=2)

plt.tight_layout()




