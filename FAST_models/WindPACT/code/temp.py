"""
compare 1.5 MW responses
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import os
import matplotlib.pyplot as plt
import numpy as np

fileID = '_00000'
turb_ID = ['C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\FAST7\\WP1500_FAST_v7','WP1500', \
        'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\FAST7\\WP1.5A08V03','WP1.5A08V03']
    
plt.figure(1)
plt.clf()    
    
for i in range(2):
    turb_dir = turb_ID[2*i]
    TName   = turb_ID[2*i+1]
    leg_str = ['FASTCert 1.5','JRink 1.5'][i]

    fname = os.path.join(turb_dir,TName+fileID+'.out')

    FAST = jr.ReadFASTFile(fname)

    Data   = FAST['Data']
    Fields = FAST['Fields']
    
    PlotFields = ['Time','WindVxi','GenSpeed','GenPwr','BldPitch2',
                   'TSR','GenTq']
    
    nplots = len(PlotFields) - 1
    for i_plot in range(nplots):
        
        t = Data[:,Fields.index('Time')]
        y = Data[:,Fields.index(PlotFields[i_plot+1])]
        
        ax = plt.subplot(nplots,1,i_plot+1)
        ax.plot(t,y,label=leg_str)
        ax.set_title(PlotFields[i_plot+1])
        ax.legend()
        
plt.tight_layout()