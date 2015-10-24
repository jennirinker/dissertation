"""
Plot steady-state response
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import os
import numpy as np
import scipy.io as scio
import matplotlib.pyplot as plt


# set directory and turbine name
#turb_dir,TName,fignum,ylim1,ylim2 = 'C:\\Users\\jrinker\\Documents\\' + \
#    'GitHub\\dissertation\\FAST_models\\FAST7\\WP0.75A08V00', \
#    'WP0.75A08V00',1,[0,2000],[-0.5,2.0]
#turb_dir,TName,fignum,ylim1,ylim2 = 'C:\\Users\\jrinker\\Documents\\' + \
#    'GitHub\\dissertation\\FAST_models\\FAST7\\WP1.5A08V03', \
#    'WP1.5A08V03',2,[0,2000],[-0.5,2.5]
#turb_dir,TName,fignum,ylim1,ylim2 = 'C:\\Users\\jrinker\\Documents\\' + \
#    'GitHub\\dissertation\\FAST_models\\FAST7\\WP1500_FAST_v7', \
#    'WP1500',3,[0,2000],[-0.5,2.5]
#turb_dir,TName,fignum,ylim1,ylim2 = 'C:\\Users\\jrinker\\Documents\\' + \
#    'GitHub\\dissertation\\FAST_models\\FAST7\\WP3.0A02V02', \
#    'WP3.0A02V02',4,[0,4000],[-0.5,2.5]
turb_dir,TName,fignum,ylim1,ylim2 = 'C:\\Users\\jrinker\\Documents\\' + \
    'GitHub\\dissertation\\FAST_models\\FAST7\\WP5.0A04V00', \
    'WP5.0A04V00',5,[0,6000],[-1.0,3.5]
    
# load steady-state dictionary
mdict = scio.loadmat(os.path.join(turb_dir,'steady_state',
                                TName+'_SS.mat'),squeeze_me=True)
LUT    = mdict['SS']
fields = [str(x).strip() for x in mdict['Fields']]

# get all values
WindVxi   = LUT[:,fields.index('WindVxi')]

# initialize figure and axes
fig1 = plt.figure(fignum,figsize=(7,9))
plt.clf()

fig1,ax1,ax2,ax3 = jr.PlotSSTurbineResponse(WindVxi,LUT,fields,fig=fig1)

ax1.set_title(TName+' Steady State')
ax1.set_ylim(ylim1)
ax3.set_ylim(ylim2)
