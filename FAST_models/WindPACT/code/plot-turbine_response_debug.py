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
fileID,fignum = '42331',1
#fileID,fignum = '93013',2


#PlotFields = ['Time','WindVxi','RotSpeed','GenPwr',
#              'BldPitch1','TSR','GenTq','TTDspFA','OoPDefl1']
PlotFields = ['Time','WindVxi','RotSpeed',
              'BldPitch1','GenTq','TTDspFA','OoPDefl1']
              
# initialize figure
fig = plt.figure(fignum,figsize=(6.5,10))
plt.clf()

# load fast file
FASTfname = TName+'_'+fileID
FASTfpath = os.path.join(turb_dir,FASTfname+'.out')
FAST = jr.ReadFASTFile(FASTfpath)

# plot results
t = FAST['Data'][:,FAST['Fields'].index('Time')]

tplot = [150,155]
mask = (t >= tplot[0]) & (t <= tplot[1])

for i_plot in range(len(PlotFields)-1):
    PlotField = PlotFields[i_plot+1]
    
    ax = fig.add_subplot(len(PlotFields)-1,1,i_plot+1)
    
    if PlotField == 'WindVxi':
        ax.plot(t[mask],11*np.ones(t[mask].shape),'k:')
    
    y = FAST['Data'][:,FAST['Fields'].index(PlotField)]
    

    ax.plot(t[mask],y[mask])
    ax.set_title(PlotFields[i_plot+1])
    
    plt.locator_params(axis='y',nbins=4)
    #        if i_plot < len(PlotFields)-2:
    #            ax.set_ylim([np.mean(y)*0.95,np.mean(y)*1.05])
#    ax.set_xlim([150,160])
#    ax.margins(y=.1)
            
plt.tight_layout()
fig.suptitle(FASTfname,x=0.01,y=0.99,fontsize='large',ha='left')