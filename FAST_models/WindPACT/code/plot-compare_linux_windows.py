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


# specify the directory to write the files to
basedir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
        'dissertation\\FAST_models\\FAST7'
tnames = ['WP0.75A08V00','WP1.5A08V03',
          'WP3.0A02V02','WP5.0A04V00']
        
fileID = '_00000'    
PlotFields = ['Time','WindVxi','RotSpeed','GenPwr',
              'BldPitch1','TSR','GenTq','TwrBsMxt','OoPDefl1']
leg_str = ['Windows','Peregrine']
c = ['b','r']   
       
for i_turb in range(len(tnames)):
    TName      = tnames[i_turb]
    TNameLinux = TName+'_Linux'

    turbs = [TName, os.path.join(basedir,TName),\
            TName, os.path.join(basedir,TNameLinux)]

    fig = plt.figure(i_turb+1,figsize=(6.5,10))
    plt.clf()
    
    for i in range(2):
        
        TName = turbs[2*i]
        turb_dir = turbs[2*i+1]
        
        fname = os.path.join(turb_dir,TName+fileID)
        FAST = jr.ReadFASTFile(fname+'.out')
        
        t = FAST['Data'][:,FAST['Fields'].index('Time')]
        
        for i_plot in range(len(PlotFields)-1):
            ax = fig.add_subplot(len(PlotFields)-1,1,i_plot+1)
            
            y = FAST['Data'][:,FAST['Fields'].index(PlotFields[i_plot+1])]
            ax.plot(t,y,c[i],label=leg_str[i])
            ax.set_title(PlotFields[i_plot+1])
    #        if i_plot < len(PlotFields)-2:
    #            ax.set_ylim([np.mean(y)*0.95,np.mean(y)*1.05])
            if i_plot == 0:
                ax.legend()
            
    plt.tight_layout()
    fig.suptitle(TName,x=0.14,y=0.99,fontsize='large')