"""
still dicking around with movie stuff
"""
"""
messing around with processed statistics
"""
import scipy.io as scio
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

# define turbine name
turb_names = ['WP0.75A08V00','WP1.5A08V03',
              'WP3.0A02V02','WP5.0A04V00']

# hard-code fields for now
fields = ['Time', 'WindVxi', 'WindVyi', 'WindVzi', 'OoPDefl1', 'IPDefl1',
          'TipDzb1', 'TwrClrnc1', 'OoPDefl2', 'IPDefl2', 'TipDzb2', 'TwrClrnc2',
          'OoPDefl3', 'IPDefl3', 'TipDzb3', 'TwrClrnc3', 'BldPitch1', 'BldPitch2',
          'BldPitch3', 'Azimuth', 'RotSpeed', 'RotAccel', 'GenSpeed', 'GenAccel',
          'TSR', 'TTDspFA', 'TTDspSS', 'TTDspAx', 
           'RootFzb1', 'RootMIP1', 'RootMOoP1', 'RootMzb1',
          'RootMEdg1', 'RootMFlp1', 'RootFzb2', 'RootMIP2', 'RootMOoP2',
          'RootMzb2', 'RootMEdg2', 'RootMFlp2', 'RootFzb3', 'RootMIP3',
          'RootMOoP3', 'RootMzb3', 'RootMEdg3', 'RootMFlp3', 'RotThrust',
          'LSSGagFya', 'LSSGagFza', 'LSSGagFys', 'LSSGagFzs', 'RotTorq',
          'CThrstArm', 'LSShftPwr', 'LSShftCq', 'LSShftCp', 'LSShftCt',
          'HSShftTq', 'HSShftPwr', 'HSShftCq', 'HSShftCp', 'GenTq', 'GenPwr',
          'GenCq', 'GenCp', 'YawBrFzp', 'YawBrMzp', 'TwrBsFxt', 'TwrBsFyt',
          'TwrBsFzt', 'TwrBsMxt', 'TwrBsMyt', 'TwrBsMzt']

# define list of parameters
us   = [5,7,9,10,10.5,11,11.5,12,13,16,19,22]
tis  = [0.1,0.2,0.3,0.4,0.5]
ls   = [10**1.5,10**2.,10**2.5,10**3]
rhos = [0.,0.1,0.2,0.3,0.4]

# define statistics and parameter to plot
#stat = 'max'
#parm = 'OoPDefl1'
#stat = 'mean'
#parm = 'GenPwr'
#stat = 'max'
#parm = 'TTDspFA'
stat = 'min'
parm = 'OoPDefl1'

# plotting options
#u_lo,u_hi = 6,11.6                          # wind velocity mask
u_lo,u_hi = 0,30                          # wind velocity mask
color_opts = ['b','r','g','m','c']          # marker colors
marker_opts = ['o','^','s','p','d']         # marker types

# intermediate parameters
n_fields = len(fields)

# define turbine
i_turb = 3

# create plot
turb_name = turb_names[i_turb]

# load data
stats_dict = scio.loadmat(turb_name+'_stats.mat',squeeze_me=True)
proc_stats = stats_dict['proc_stats']
calc_stats = [s.rstrip() for s in stats_dict['calc_stats']]
fnames     = [s.rstrip() for s in stats_dict['fnames']]
#fields     = [s.rstrip() for s in stats_dict['fields']]

# calculate specified statistic
zs   = proc_stats[:,calc_stats.index(stat)*n_fields + \
                                    fields.index(parm)]

# extract corresponding mean wind speed and ti
xs = np.empty(zs.size)
ys = np.empty(zs.size)
rhos_idx = np.empty(zs.size)
for i_f in range(zs.size):
    file_id =  fnames[i_f].rstrip('.out').split('_')[1]
    xs[i_f]  = us[int(file_id[0],16)]
    ys[i_f]   = tis[int(file_id[1],16)]
    rhos_idx[i_f] = int(file_id[3],16)

# mask data to only plot certain portion
x_mask = np.ma.masked_outside(xs,u_lo,u_hi)
y_mask = np.ma.masked_where(np.ma.getmask(x_mask),ys)
z_mask = np.ma.masked_where(np.ma.getmask(x_mask),zs)
rhos_mask = np.ma.masked_where(np.ma.getmask(x_mask),rhos_idx)

# initialize figure
fig = plt.figure(i_turb+1,figsize=(6,6))
plt.clf()
ax = fig.add_subplot(111, projection='3d')

# create scatterplot
for i_rho in range(5):
    idx = np.where(rhos_mask == i_rho)
    xplot = x_mask[idx]
    yplot = y_mask[idx]
    zplot = z_mask[idx]
    
    ax.scatter(xplot, yplot, zplot, s=9,
               c=color_opts[i_rho], 
               marker=marker_opts[i_rho],
               edgecolors='none',
               label=r'$\rho_{:d}$'.format(i_rho))    
                     
# prettify
ax.set_xlabel('Mean wind speed')
ax.set_ylabel('Turbulence intensity')
ax.set_zlabel('{:s} of {:s}'.format(stat,parm))
fig.suptitle('{:s}: {:s} of {:s}'.format(turb_name,stat,parm),fontsize='large')
plt.tight_layout()

# rotate and save pngs
angles = np.arange(0,360,4)
elev = 15

for i in range(angles.size):
    ax.view_init(elev, angles[i])
    plt.draw()
    fig.savefig('pics\\{:s}_{:d}.png'.format(turb_name,i))
    
print('done')