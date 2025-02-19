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
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP0.75A08V00','WP0.75A08V00'
turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\FAST7\\WP0.75A08V00_newGBR','WP0.75A08V00'
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

out_fnames = [f for f in os.listdir(turb_dir) if f.endswith('.out')]

for i_f in range(len(out_fnames)):
    
    # initialize figure
    fig1 = plt.figure(i_f+1,figsize=(6.5,10))
    plt.clf()

    # load fast file
    FASTfname = out_fnames[i_f]
    FASTfpath = os.path.join(turb_dir,FASTfname)
    FAST = jr.ReadFASTFile(FASTfpath)
    
    # plot results
    t = FAST['Data'][:,FAST['Fields'].index('Time')]
    jr.PlotTurbineResponse(t,FAST['Data'],FAST['Fields'],fig=fig1)
    
    fig1.suptitle(FASTfname,fontsize='large')