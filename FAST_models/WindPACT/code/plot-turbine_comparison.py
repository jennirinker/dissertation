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
#turbs = ['WP1500', os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#        'dissertation\\FAST_models\\FAST7','WP1500_FAST_v7'),\
#        'WP1.5A08V03', os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#        'dissertation\\FAST_models\\FAST7','WP1.5A08V03')]
#turbs = ['WP1.5A08V03', os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#        'dissertation\\FAST_models\\FAST7','WP1.5_Linux'),\
#        'WP1.5A08V03', os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#        'dissertation\\FAST_models\\FAST7','WP1.5A08V03')]
#turbs = ['WP0.75A08V00', os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#        'dissertation\\FAST_models\\FAST7','WP0.75A08V00'),\
#        'WP1.5A08V03', os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#        'dissertation\\FAST_models\\FAST7','WP1.5A08V03')]
#turbs = ['WP0.75A08V00', os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#        'dissertation\\FAST_models\\FAST7','WP0.75A08V00'),\
#        'WP0.75A08V00', os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#        'dissertation\\FAST_models\\FAST7','WP0.75A08V00')]
turbs = ['WP0.75A08V00', os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
        'dissertation\\FAST_models\\FAST7','WP0.75A08V00'),\
        'WP0.75A08V00', os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
        'dissertation\\FAST_models\\FAST7','WP0.75A08V00_dynin')]
#turbs = ['WP0.75A08V00', os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#        'dissertation\\FAST_models\\FAST7','WP0.75A08V00'),\
#        'WP0.75A08V00', os.path.join('\\\\monsoon-data\\Public\\JRinker\\' + \
#        'fast_simulations\\FastDir\\SmallRun','WP0.75A08V00')]
#        'WP0.75A08V00', os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#        'dissertation\\FAST_models\\FAST7','WP0.75A08V00_stifftwr'),\
#        'WP0.75A08V00', os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#        'dissertation\\FAST_models\\FAST7','WP0.75A08V00_stiffblds')]
#fileIDs = ['_42331','_42331','_42331','_42331','_42331']
#fileIDs = ['_91242','_91242','_91242','_91242','_91242']
fileIDs = ['_00000','_00000','_00000','_00000','_00000']
#fileIDs = ['_00000','_00001']

t_plot = [30,630]

PlotFields = ['Time','WindVxi','RotSpeed','GenPwr',
              'BldPitch1','TSR','GenTq','YawBrTAxp','OoPDefl1']
#leg_str = ['Linux 1.5','Windows 1.5']
#leg_str = ['Standard pitch','Modified GS']
#leg_str = ['Standard model','Modified GBR','Stiffened tower',
#           'Stiffened blades','Extra damp']
leg_str = ['EQUIL','DYNIN']
#leg_str = ['Desktop','Monsoon']
c = ['b','r','g','c','m','y','k']

fig1 = plt.figure(5,figsize=(6.5,10))
plt.clf()

for i in range(len(turbs)/2):
    
    TName = turbs[2*i]
    turb_dir = turbs[2*i+1]
    
    fname = os.path.join(turb_dir,TName+fileIDs[i])
    FAST = jr.ReadFASTFile(fname+'.out')
    
    t = FAST['Data'][:,FAST['Fields'].index('Time')]
    
    print(FAST['Data'][:,FAST['Fields'].index('OoPDefl1')].max())
    
    for i_plot in range(len(PlotFields)-1):
        ax = fig1.add_subplot(len(PlotFields)-1,1,i_plot+1)
        
        y = FAST['Data'][:,FAST['Fields'].index(PlotFields[i_plot+1])]
        ax.plot(t,y,c[i],label=leg_str[i])
        ax.set_title(PlotFields[i_plot+1])
#        if i_plot < len(PlotFields)-2:
#            ax.set_ylim([np.mean(y)*0.95,np.mean(y)*1.05])
        if i_plot == 0:
            ax.legend()
            
        ax.set_xlim(t_plot)
        
plt.tight_layout()