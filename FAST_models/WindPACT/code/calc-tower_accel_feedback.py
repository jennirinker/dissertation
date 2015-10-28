"""
investigating tower acceleration feedback
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
import os
import matplotlib.pyplot as plt


# set directory and turbine name
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub' + \
#        '\\dissertation\\FAST_models\\verification\\WP0.75A08V00','WP0.75A08V00'
turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\FAST7\\WP0.75A08V00','WP0.75A08V00'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP1.5A08V03','WP1.5A08V03'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP3.0A02V02','WP3.0A02V02'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP5.0A04V00','WP5.0A04V00'
fileID,fignum = '42331',11
#fileID,fignum = '93013',2


#PlotFields = ['Time','WindVxi','RotSpeed','GenPwr',
#              'BldPitch1','TSR','GenTq','TTDspFA','OoPDefl1']
PlotFields = ['Time','WindVxi','RotSpeed',
              'BldPitch1','GenTq','YawBrTAxp','OoPDefl1']
              
# initialize figure
fig = plt.figure(fignum,figsize=(6.5,10))
plt.clf()

# load fast file
FASTfname = TName+'_'+fileID
FASTfpath = os.path.join(turb_dir,FASTfname+'.out')
FAST = jr.ReadFASTFile(FASTfpath)

# plot results
t = FAST['Data'][:,FAST['Fields'].index('Time')]

tplot = [150,160]
mask = (t >= tplot[0]) & (t <= tplot[1])

# plot wind

PlotField = 'YawBrTAxp'
ax = fig.add_subplot(3,1,1)
y = FAST['Data'][:,FAST['Fields'].index(PlotField)]
ax.plot(t[mask],y[mask])
ax.set_title(PlotField)
plt.locator_params(axis='y',nbins=4)
ax.grid('on')

PlotField = 'BldPitch1'
#ax2 = ax.twinx()
ax = fig.add_subplot(3,1,2)
y = FAST['Data'][:,FAST['Fields'].index(PlotField)]
ax.plot(t[mask],y[mask],'r')
ax.set_ylim([0,15])

ax = fig.add_subplot(3,1,3)
ax.plot(t[mask],y[mask]+0.3*FAST['Data'][:,FAST['Fields'].index('YawBrTAxp')][mask])
ax.set_ylim([0,15])
