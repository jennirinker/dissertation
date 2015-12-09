"""
Summarize processed wind parameters from a given dataset
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import os
import numpy as np
import matplotlib.pyplot as plt

# style file
plt.style.use('C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
            '2016-02-15_dissertation\\figure_code\\duke_paper.mplstyle')

# choose which dataset
i_dat     = 2
datasets  = ['NREL','fluela','CM06']
dataset   = datasets[i_dat]

# define directory where wind parameters are stored (unused for matlab)
basedir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                  'processed_data'
                  
print('\nProcessing dataset {}'.format(dataset))

if ('mat' in dataset):
    fields, raw_parms = jr.loadNRELmatlab()
    dataset = 'NREL'
else:
    fname = dataset + '-metadata.mat'
    fpath = os.path.join(basedir,fname)
    fields, raw_parms = jr.loadmetadata(fpath)

# get fieldnames for instrument ID and wind direction
fieldname_ID  = [s for s in fields if ('Height' in s or 'ID' in s)][0]
fieldname_dir = [s for s in fields if ('Direction' in s)][0]
fieldname_spd = [s for s in fields if ('Cup' in s)][0]

# custom figuze sizes for different datasets
if (dataset == 'NREL'):
    fignum = 1
    figsize = (6.5,5.)
    nax_y,nax_x = 2, 3
elif (dataset == 'fluela'):
    fignum = 2
    figsize = (6.5,2.5)
    nax_y,nax_x = 1,3
elif (dataset == 'CM06'):
    fignum = 3
    figsize = (6.5,6.0)
    nax_y,nax_x = 3,4
    print('DEFAULTING TO ROTATING MANUALLY FOR VIS.')
else:
    KeyError('Dataset not coded')
    
# instrument names
IDs = jr.datasetSpecs(dataset)[2]
n_instr = len(IDs)

# inititalize figure
plt.figure(fignum,figsize=figsize)
plt.clf()
plt.suptitle('Windroses for Dataset {:s}'.format(dataset))

# plotting properties
plot_perc = 0.95
wd_perc, ht_perc = 0.75, 0.75

dax_x,dax_y = plot_perc/nax_x, plot_perc/nax_y
wax,  hax   = dax_x*wd_perc,dax_y*ht_perc

offst_x = 1-plot_perc + wax*(1-wd_perc)/2.
offst_y = (1-plot_perc)/2.
axs_x = np.arange(offst_x,1.,dax_x)
axs_y = np.arange(offst_y,1.,dax_y)

nbins_theta   = 30
nbins_windspd = 20
bottom = 0
width = 0.7 * (2*np.pi) / nbins_theta

print('**\nNEED TO ADD COLORBAR\n')

# loop through instruments
for i_id in range(n_instr):
    
    # isolate data from that instrument
    ID     = IDs[i_id]
    data   = raw_parms[raw_parms[:,fields.index(fieldname_ID)] == ID,:]
    n_data = data.shape[0]
    
    
    # get wind spped and direction histogram info
    wind_spd = data[:,fields.index(fieldname_spd)]
    wind_dir = data[:,fields.index(fieldname_dir)] % 360
    if (dataset == 'CM06'):
        wind_dir = (120 - wind_dir)%360
        
    # remove NaN values
    idcs_notnan = np.logical_and(np.logical_not(np.isnan(wind_spd)),
                                 np.logical_not(np.isnan(wind_dir)))
    wind_spd = wind_spd[idcs_notnan]
    wind_dir = wind_dir[idcs_notnan]
    
#    radii,bin_edges = np.histogram(wind_dir,
#                                   bins=nbins)
#    theta           = (0.5*bin_edges[1:] + 0.5*bin_edges[:-1])*np.pi/180
    counts,bintheta_edges,binspd_edges = np.histogram2d(wind_dir,wind_spd,
                                   bins=[nbins_theta,nbins_windspd])
    theta = (0.5*bintheta_edges[1:] + 0.5*bintheta_edges[:-1])*np.pi/180
    spds  = (0.5*binspd_edges[1:] + 0.5*binspd_edges[:-1])
    
    # initialize
    iax_x = i_id % nax_x
    iax_y = i_id / nax_x
    ax = plt.axes([axs_x[iax_x],axs_y[iax_y],wax,hax], polar=True)

    bottom = np.zeros(nbins_theta)
    for i_spd in range(len(binspd_edges)-1):
        spd       = spds[i_spd]
        count_spd = counts[:,i_spd]
        bars = ax.bar(theta, count_spd, width=width, bottom=bottom)
        bottom += count_spd
    
        # prettify axes
        for r, bar in zip(count_spd, bars):
            bar.set_facecolor(plt.cm.Reds(spd/float(spds.max())))
            bar.set_edgecolor(plt.cm.Reds(spd/float(spds.max())))
#            bar.set_edgecolor('none')
        ax.set_theta_offset(0.5*np.pi)
        ax.set_theta_direction(-1)
        ax.set_yticklabels([])
    
    # add title to axes
    if (dataset in ('NREL','fluela')):
        ax_title = 'Height: {:d} m'.format(ID)
    elif (dataset in ('CM06')):
        ax_title = 'Sonic ID: {:d}'.format(ID)
    ax.set_title(ax_title,y=1.12)

    
plt.show()

