"""
messing around with processed statistics
"""
import scipy.io as scio
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm

# define turbine name
turb_names = ['WP0.75A08V00','WP1.5A08V03',
              'WP3.0A02V02','WP5.0A04V00']

# define list of parameters
us   = [5,7,9,10,10.5,11,11.5,12,13,16,19,22]
tis  = [0.1,0.2,0.3,0.4,0.5]
ls   = [10**1.5,10**2.,10**2.5,10**3]
rhos = [0.,0.1,0.2,0.3,0.4]

U_Ti = 1    # 0: L/rho, 1: U/Ti

# define statistics and parameter to plot
stat = 'min'
parm = 'OoPDefl1'
#stat = 'max'
#parm = 'RootMFlp1'
#stat = 'max'
#parm = 'TTDspFA'
#stat = 'max'
#parm = 'OoPDefl1'

# plotting options
savefig = 0
#u_lo,u_hi = 8.5,10.5                          # wind velocity mask
if U_Ti:
#    x_lo,x_hi = 0,30                          # wind velocity mask
    x_lo,x_hi = 0,9                          # wind velocity mask
else:
    x_lo, x_hi = 1, 4
cmap = cm.Reds
elev, az = 1,-90

# loop through turbines
for i_turb in range(len(turb_names)):
    turb_name = turb_names[i_turb]

    # load data
    stats_dict = scio.loadmat(turb_name+'_stats.mat',squeeze_me=True)
    proc_stats = stats_dict['proc_stats']
    calc_stats = [s.rstrip() for s in stats_dict['calc_stats']]
    fnames     = [s.rstrip() for s in stats_dict['fnames']]
    fields     = [s.rstrip() for s in stats_dict['fields']]
    n_fields = len(fields)

    # calculate specified statistic
    zs   = proc_stats[:,calc_stats.index(stat)*n_fields + \
                                        fields.index(parm)]
    
#    # print turbine index and filename for max stat
#    max_idcs = zs.argsort()[-5:][::-1]
#    print(i_turb)
#    print(zs[max_idcs])
#    print([fnames[i] for i in max_idcs])
    
    # extract corresponding mean wind speed and ti
    xs    = np.empty(zs.size)
    ys    = np.empty(zs.size)
    c_idx = np.empty(zs.size)
    for i_f in range(zs.size):
        file_id =  fnames[i_f].rstrip('.out').split('_')[1]
        if U_Ti:
            xs[i_f]  = us[int(file_id[0],16)]
            ys[i_f]   = tis[int(file_id[1],16)]
            c_idx[i_f] = int(file_id[3],16)
            n_c = len(rhos)
            xlabel = 'Mean wind speed'
            ylabel = 'Turbulence intensity'
        else:
            xs[i_f]  = np.log10(ls[int(file_id[2],16)])
            ys[i_f]   = rhos[int(file_id[3],16)]
            c_idx[i_f] = int(file_id[0],16)
            n_c = len(us)
            xlabel = 'log$_{10}$(L)'
            ylabel = 'Rho'

    # mask data to only plot certain portion
    x_mask = np.ma.masked_outside(xs,x_lo,x_hi)
    y_mask = np.ma.masked_where(np.ma.getmask(x_mask),ys)
    z_mask = np.ma.masked_where(np.ma.getmask(x_mask),zs)
    c_mask = np.ma.masked_where(np.ma.getmask(x_mask),c_idx)

    # initialize figure
    fig = plt.figure(i_turb+1,figsize=(6,6))
    plt.clf()
    ax = fig.add_subplot(111, projection='3d')

    # create scatterplot
    for i_c in range(n_c):
        idx = np.where(c_mask == i_c)
        xplot = x_mask[idx]
        yplot = y_mask[idx]
        zplot = z_mask[idx]
        
#        ax.scatter(xplot, yplot, zplot, s=9,
#                   c=cmap(i_c/float(n_c)), 
#                   edgecolors='none',
#                   label=r'$\rho_{:d}$'.format(i_c))  
        ax.scatter(xplot, yplot, zplot, s=9,
                   c='k', 
                   edgecolors='none',
                   label=r'$\rho_{:d}$'.format(i_c))    
                      
    # prettify
    ax.view_init(elev,az)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel('{:s} of {:s}'.format(stat,parm))
    fig.suptitle('{:s}: {:s} of {:s}'.format(turb_name,stat,parm),fontsize='large')
    plt.tight_layout()
        
    # save handle for funzies
    if savefig:
        fig.savefig('C:\\Users\\jrinker\\Desktop\\resp_{:d}.png'.format(i_turb))
    
    