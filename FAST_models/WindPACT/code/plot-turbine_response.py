"""
Plot comparison of turbine models
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
import os
import matplotlib.pyplot as plt


# set directory and turbine name
turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\FAST7\\WP0.75A08V00','WP0.75A08V00'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP0.75A08V00_newGBR','WP0.75A08V00'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP0.75A08V00_newGBR','WP0.75A08V00'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP0.75A08V00_steady','WP0.75A08V00'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP1500_FAST_v7','WP1500'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP1.5A08V03','WP1.5A08V03'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP3.0A02V02','WP3.0A02V02'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP5.0A04V00','WP5.0A04V00'
#fileID = '00000'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub' + \
#        '\\dissertation\\FAST_models\\verification\\WP0.75A08V00','WP0.75A08V00'
fileID = ['24134', '24142', '42331', '91242'][0]
#fileID = ['24211', '24214', '24220'][0]
#fileID = ['24322', '24324', '24343'][0]
t_lim = [30,300]       # WP0.75
#t_lim = [240,300]      # WP1.5
grid_flg = 'on'

# define turbine version
ext,fignum = '',1
#ext,fignum = '_newGBR',2
#ext,fignum = '_stiffblds',3
#ext,fignum = '_equil',4
turb_dir += ext

# initialize figure
fig1 = plt.figure(fignum,figsize=(6.5,10))
plt.clf()

# load fast file
FASTfname = TName+'_'+fileID
FASTfpath = os.path.join(turb_dir,FASTfname+'.out')
FAST = jr.ReadFASTFile(FASTfpath)

# plot results
t = FAST['Data'][:,FAST['Fields'].index('Time')]
fig,ax1,ax2,ax3 = jr.PlotTurbineResponse(t,FAST['Data'],FAST['Fields'],fig=fig1)

#ax3.set_ylim([-10,6])

ax1.set_xlim(t_lim)
ax2.set_xlim(t_lim)
ax3.set_xlim(t_lim)
ax1.grid(grid_flg)
ax2.grid(grid_flg)
ax3.grid(grid_flg)

fig1.suptitle(FASTfname+ext,fontsize='large')
