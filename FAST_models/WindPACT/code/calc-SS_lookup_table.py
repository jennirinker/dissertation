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

# initialize look-up table
fSSname = TName + '_SS.mat'
fSSpath = os.path.join(turb_dir,'steady_state',fSSname)
saveFields = ['WindVxi',
          'GenSpeed','LSShftPwr','GenPwr','RotThrust','RotTorq',
          'RotSpeed','BldPitch1','GenTq','TSR',
          'OoPDefl1','IPDefl1','TTDspFA','TTDspSS']
LUT = np.empty((0,len(saveFields)))

# initial parameters
TMax      = 120.0
T_ss      = 60.0

# change directory to turbine directory to run FAST
os.chdir(turb_dir)

# set first initial conditions
BlPitch0  = 2.6*np.ones(3)
RotSpeed0 = 6.0

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
               BlPitch0=BlPitch0,RotSpeed0=RotSpeed0,
               wind_dir=wind_dir,fileID=fileID,TMax=TMax)
               
    # run FAST
    print('Processing wind speed {:.1f}'.format(wind_speed))
    FASTfname = TName+'_'+fileID
    os.system('FAST.exe '+FASTfname+'.fst')
                  
    # load FAST files
    FAST = jr.ReadFASTFile(FASTfname+'.out')
    
    # loop through and save steady-state values
    n_t = FAST['Data'][:,FAST['Fields'].index('Time')].size*T_ss/TMax
    row = np.empty(len(saveFields))
    for i_parm in range(len(saveFields)):
        parm = saveFields[i_parm]
        x = FAST['Data'][:,FAST['Fields'].index(parm)]
                      
        # calculate and save last value
        x_SS = np.mean(x[-n_t:])
        row[i_parm] = np.array(x_SS)
    LUT = np.vstack((LUT,row.reshape((1,row.size))))
    
    # rearrange to increasing wind speed
    LUT = LUT[LUT[:,saveFields.index('WindVxi')].argsort()]
                  
    # set initial conditions for next round
    BlPitch0  = LUT[-1,saveFields.index('BldPitch1')]*np.ones(3)
    RotSpeed0 = LUT[-1,saveFields.index('RotSpeed')]
    
    # plot results and save figure
    t = FAST['Data'][:,FAST['Fields'].index('Time')]
    jr.PlotTurbineResponse(t,FAST['Data'],FAST['Fields'],fig=fig1)
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
LUT = LUT[LUT[:,0].argsort()]
   
# save LUT   
mdict = {}
mdict['SS'] = LUT   
mdict['Fields'] = saveFields 
scio.savemat(fSSpath,mdict)
print('Look-up table saved.')