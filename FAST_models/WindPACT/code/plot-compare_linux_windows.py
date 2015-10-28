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
#basedir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#        'dissertation\\FAST_models\\FAST7'
tnames = ['WP0.75A08V00','WP1.5A08V03',
          'WP3.0A02V02','WP5.0A04V00']
basedir = 'C:\\Users\\jrinker\\Documents\\GitHub' + \
        '\\dissertation\\FAST_models\\verification'
#tnames = ['WP0.75A08V00']
        
fileID = '_42331'    
PlotFields = ['Time','WindVxi','RotSpeed','GenPwr',
              'BldPitch1','TSR','GenTq','TwrBsMxt','OoPDefl1']
leg_str = ['Windows','Peregrine']
c = ['b','r']   
       
# loop through the turbines
fig_idx = 1
for i_turb in range(len(tnames)):
    TName      = tnames[i_turb]
#    TNameLinux = TName+'_Linux'
    TNameLinux = TName+'_Peregrine'
    
    turbs = [TName, os.path.join(basedir,TName),\
            TName, os.path.join(basedir,TNameLinux)]
            
    # loop through the .fst files
    fst_names = [fname for fname in os.listdir(turbs[1]) \
                if fname.endswith('.out')]
    for i_fst in range(len(fst_names)):
        
        fst_name = fst_names[i_fst]

        fig = plt.figure(fig_idx,figsize=(6.5,10))
        plt.clf()
    
        # for windows and peregrine
        for i in range(2):
            
            TName = turbs[2*i]
            turb_dir = turbs[2*i+1]
            
            fst_path = os.path.join(turb_dir,fst_name)
            FAST = jr.ReadFASTFile(fst_path)
            
            t = FAST['Data'][:,FAST['Fields'].index('Time')]
            
            for i_plot in range(len(PlotFields)-1):
                PlotField = PlotFields[i_plot+1]
                
                ax = fig.add_subplot(len(PlotFields)-1,1,i_plot+1)
                
                if PlotField == 'WindVxi':
                    ax.plot(t,11*np.ones(t.shape),'k:')
                
                y = FAST['Data'][:,FAST['Fields'].index(PlotField)]
                ax.plot(t,y,c[i],label=leg_str[i])
                ax.set_title(PlotFields[i_plot+1])
        #        if i_plot < len(PlotFields)-2:
        #            ax.set_ylim([np.mean(y)*0.95,np.mean(y)*1.05])
                if i_plot == 0:
                    ax.legend()
                    
                ax.set_xlim([30,630])
                
        plt.tight_layout()
        fig.suptitle(fst_name,x=0.01,y=0.99,fontsize='large',ha='left')
#        plt.savefig('C:\\Users\\jrinker\\Desktop\\'+fst_name[:-4]+'.png')
        
        fig_idx += 1