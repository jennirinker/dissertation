"""
Get steady-state look-up table for a turbine
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import os
import numpy as np
import scipy.io as scio
import matplotlib.pyplot as plt

# wind speeds to calculate steady-state
# *************** MUST BE INCREASING ORDER ***************
wind_speeds = np.arange(3,25,0.25)

# directory where wind files are located
wind_dir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\wind_files'

# set directory and turbine name
turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\FAST7\\WP0.75A08V00','WP0.75A08V00'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP1500_FAST_v7','WP1500'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP1.5A08V03','WP1.5A08V03'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP3.0A02V02','WP3.0A02V02'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP5.0A04V00','WP5.0A04V00'

# initialize look-up table
fSSname = TName + '_SS.mat'
fSSpath = os.path.join(turb_dir,'steady_state',fSSname)

# initial parameters
TMax      = 140.0
T_ss      = 80.0

# change directory to turbine directory to run FAST
os.chdir(turb_dir)

# set first initial conditions
ICDict = {}
ICDict['BlPitch(1)'],ICDict['BlPitch(2)'],ICDict['BlPitch(3)'] = 2.6,2.6,2.6
ICDict['OoPDefl'],ICDict['IPDefl'],ICDict['TeetDefl'] = 0.,0.,0.
ICDict['Azimuth'],ICDict['RotSpeed'],ICDict['NacYaw'] = 0.,6.,0.
ICDict['TTDspFA'],ICDict['TTDspSS'] = 0.,0.
ICDict['TMax'],ICDict['TStart'] = TMax, 0.

# initialize stuff for plotting
PlotFields = ['WindVxi','RotSpeed','GenPwr',
              'BldPitch1','TSR','GenTq','TwrBsMxt'] 
fig1 = plt.figure(1,figsize=(6.5,10))
plt.clf()

# loop through wind speeds
for i_WS in range(wind_speeds.size):
    wind_speed = wind_speeds[i_WS]
        
    # set fileID for FAST run
    fileID     = '{:.0f}'.format(i_WS).zfill(5)
    
    # create wind filename
    wind_fname = 'NoShr_'+'{:2.1f}'.format(wind_speed).zfill(4)+'.wnd'
    
    # check if wind file exists, make it if not
    wind_fpath = os.path.join(wind_dir,wind_fname)
    if not os.path.exists(wind_fpath):
        jr.writeSteadyWind(wind_speed,wind_dir=wind_dir)
            
    # create FAST files
    jr.writeFASTFiles(turb_dir,TName,wind_fname,
               wind_dir=wind_dir,fileID=fileID,
               **ICDict)
               
    # run FAST
    print('Processing wind speed {:.1f}'.format(wind_speed))
    FASTfname = TName+'_'+fileID
    os.system('FAST.exe '+FASTfname+'.fst')
                  
    # load FAST files
    FAST   = jr.ReadFASTFile(FASTfname+'.out')
    Fields = FAST['Fields']
    Data   = FAST['Data']
    
    # initialize LUT if it doesn't exist
    if i_WS == 0: LUT = np.empty((len(wind_speeds),len(Fields)))
    
    # loop through and save steady-state values
    n_t = Data[:,Fields.index('Time')].size*T_ss/TMax
    for i_parm in range(len(Fields)):
        
        # get data
        parm = Fields[i_parm]
        x = Data[:,Fields.index(parm)]
                      
        # calculate and save last value
        x_SS = np.mean(x[-n_t:])
        LUT[i_WS,i_parm] = x_SS
                              
    # set initial conditions for next round
    for key in ICDict.keys():
        if key in Fields:
            ICDict[key] = LUT[i_WS,Fields.index(key)]
        elif key+'1' in Fields:
            ICDict[key] = LUT[i_WS,Fields.index(key+'1')]
    
    # plot results and save figure
    t = Data[:,Fields.index('Time')]
    jr.PlotTurbineResponse(t,Data,Fields,fig=fig1)
    figpath = os.path.join(turb_dir,'steady_state',
                           '{:2.1f}'.format(wind_speed).zfill(4)+'.png')
    fig1.savefig(figpath)
    fig1.clf()
    
    # delete input files, move to steady-state
    os.system('del '+FASTfname+'.fst')
    os.system('del '+FASTfname+'.out')
    os.system('del '+FASTfname+'_AD.ipt')

 
plt.close(fig1)  
print('Simulations completed.')    
   
# rearrange to increasing wind speed
LUT = LUT[LUT[:,Fields.index('WindVxi')].argsort()]
   
# save LUT   
mdict = {}
mdict['SS'] = LUT   
mdict['Fields'] = Fields 
scio.savemat(fSSpath,mdict)
print('Look-up table saved.')