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
fileID,fignum = '42331',1
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

tplot = [0,160]
mask = (t >= tplot[0]) & (t <= tplot[1])

# plot wind

PlotField = 'WindVxi'
ax = fig.add_subplot(5,1,1)
ax.plot(t[mask],11*np.ones(t[mask].shape),'k--')
y = FAST['Data'][:,FAST['Fields'].index(PlotField)]
ax.plot(t[mask],y[mask],'b')
ax.set_title(PlotField)
plt.locator_params(axis='y',nbins=4)
ax.grid('on')

PlotField = 'GenTq'
ax = fig.add_subplot(5,1,2)
y = FAST['Data'][:,FAST['Fields'].index(PlotField)]
ax.plot(t[mask],y[mask])
ax.set_title(PlotField)
plt.locator_params(axis='y',nbins=4)
ax.grid('on')

PlotField = 'RotSpeed'
ax = fig.add_subplot(5,1,3)
ax.plot(t[mask],28.648*np.ones(t[mask].shape),'k--')
y = FAST['Data'][:,FAST['Fields'].index(PlotField)]
ax.plot(t[mask],y[mask])
ax.set_title('RotSpeed')
plt.locator_params(axis='y',nbins=4)
for tl in ax.get_yticklabels():
    tl.set_color('b')
ax.set_ylim([24,30])
ax.grid('on')

PlotField = 'YawBrTAxp'
ax = fig.add_subplot(5,1,4)
y = FAST['Data'][:,FAST['Fields'].index(PlotField)]
ax.plot(t[mask],y[mask])
ax.set_title(PlotField)
plt.locator_params(axis='y',nbins=4)
ax.grid('on')

PlotField = 'BldPitch1'
ax2 = ax.twinx()
y = FAST['Data'][:,FAST['Fields'].index(PlotField)]
ax2.plot(t[mask],y[mask],'r')
ax2.locator_params(axis='y',nbins=4)
for tl in ax2.get_yticklabels():
    tl.set_color('r')
ax2.set_ylim([0,15])
#ax2.plot(t[mask],y[mask]+0.3*FAST['Data'][:,FAST['Fields'].index('YawBrTAxp')][mask])

PlotField = 'OoPDefl1'
ax = fig.add_subplot(5,1,5)
y = FAST['Data'][:,FAST['Fields'].index(PlotField)]
ax.plot(t[mask],y[mask])
ax.set_title(PlotField)
plt.locator_params(axis='y',nbins=4)
ax.grid('on')

            
plt.tight_layout()
fig.suptitle(FASTfname,x=0.01,y=0.99,fontsize='large',ha='left')