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
import matplotlib as mpl

# style file
plt.style.use('C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
            '2016-02-15_dissertation\\figure_code\\duke_paper.mplstyle')

# choose which dataset
i_dat     = 0
datasets  = ['NREL','fluela','PM06']
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
elif (dataset == 'PM06'):
    fignum = 3
    figsize = (6.5,6.0)
    nax_y,nax_x = 3,4
else:
    KeyError('Dataset not coded')
    
# instrument names
IDs = jr.datasetSpecs(dataset)['IDs']
n_instr = len(IDs)

# inititalize figure
plt.figure(fignum,figsize=figsize)
plt.clf()

# plotting properties
plot_perc_x,plot_perc_y = 0.97, 1 - 0.6/figsize[1]
wd_perc, ht_perc = 0.75, 0.75

dax_x,dax_y = plot_perc_x/nax_x, plot_perc_y/nax_y
wax,  hax   = dax_x*wd_perc,dax_y*ht_perc

offst_x = 1-plot_perc_x + wax*(1-wd_perc)/2.
offst_y = (1-plot_perc_y)/4.
axs_x = np.arange(offst_x,1.,dax_x)
axs_y = np.arange(offst_y,1.,dax_y)[-2::-1]

nbins_theta   = 30
nbins_windspd = 20
bottom = 0
width = 0.7 * (2*np.pi) / nbins_theta

cb_ht,cb_dy = 0.10,-0.20
cb_ax = [0.5,1+cb_dy/figsize[1],0.4,cb_ht/figsize[1]]
cmap  = plt.cm.Reds

# get theta, wind speed, ID data
ID_data   = raw_parms[:,fields.index(fieldname_ID)]
wind_spds = raw_parms[:,fields.index(fieldname_spd)]
wind_dirs = raw_parms[:,fields.index(fieldname_dir)] % 360

# remove potential nan values
idcs_notnan = np.logical_and(np.logical_not(np.isnan(wind_spds)),
                             np.logical_not(np.isnan(wind_dirs)))
ID_data   = ID_data[idcs_notnan]
wind_spds = wind_spds[idcs_notnan]
wind_dirs = wind_dirs[idcs_notnan]

# get maximal wind speed for plotting
ws_max = np.ceil(np.nanmax(wind_spds))

# loop through instruments
for i_id in range(n_instr):
    
    # get wind speed and direction for that instrument
    ID     = IDs[i_id]
    wind_spd = wind_spds[ID_data == ID]
    wind_dir = wind_dirs[ID_data == ID]
    if (dataset == 'PM06'):
        print('DEFAULTING TO ROTATING MANUALLY FOR VIS.')
        wind_dir = (120 - wind_dir) % 360
        

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
    
        # color bars by wind speed
        for r, bar in zip(count_spd, bars):
            bar.set_facecolor(cmap(spd/float(ws_max)))
            bar.set_edgecolor(cmap(spd/float(ws_max)))
    
    # prettify axes
    ax.set_theta_offset(0.5*np.pi)
    ax.set_theta_direction(-1)
    ax.set_yticklabels([])
    if (dataset in ('NREL','fluela')):
        ax_title = 'Height: {:d} m'.format(ID)
    elif (dataset in ('PM06')):
        ax_title = 'Sonic ID: {:d}'.format(ID)
    ax.set_title(ax_title,y=1.12)

    
# add figure title and colorbar    
plt.suptitle('Wind roses for Dataset {:s}'.format(dataset),x=0.25)
axcb = plt.axes(cb_ax)
bounds = np.linspace(0,ws_max,nbins_windspd+1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
#norm = mpl.colors.Normalize(vmin=0, vmax=ws_max)
#cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap,
#                                norm=norm,
#                                orientation='horizontal')
axcb.text(1.03, 0.5,'m s$^{-1}$',
          fontsize='small',
          verticalalignment='center',
           transform=axcb.transAxes)
    
# show figure
plt.show()

