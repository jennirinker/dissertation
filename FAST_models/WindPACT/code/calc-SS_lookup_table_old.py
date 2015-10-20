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
import warnings

# wind speeds to calculate steady-state
wind_speeds = np.arange(4,25,2)

# directory where wind files are located
wind_dir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\wind_files'

# set directory and turbine name
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP1500_FAST_v7','WP1500'
turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\FAST7\\WP1.5A08V03','WP1.5A08V03'

# load look-up table if available
fLUTname = TName + '_ICs.mat'
fLUTpath = os.path.join(turb_dir,fLUTname)
if os.path.exists(fLUTpath):
    LUT = scio.loadmat(fLUTpath,squeeze_me=True)['LUT']
else:
    LUT = np.empty((0,3))

# specify whether to override values in look-up table if they already exist
overwrite = 1

# initial parameters
TMax      = 120.0
BlPitch0  = 2.6*np.ones(3)
RotSpeed0 = 10.0

os.chdir(turb_dir)

# loop through wind speeds
for i_WS in range(wind_speeds.size):
    
    wind_speed = wind_speeds[i_WS]
    
    if ((not np.any(LUT[:,0] == wind_speed)) or (overwrite)):
    
        fileID     = '{:.0f}'.format(i_WS).zfill(5)
        
        # create wind filename
        wind_fname = 'NoShr_'+'{:2.1f}'.format(wind_speed).zfill(4)+'.wnd'
        
        # check if wind file exists, make it if not
        wind_fpath = os.path.join(wind_dir,wind_fname)
        if not os.path.exists(wind_fpath):
            with open(wind_fpath,'w') as f:
                f.write('! Wind file for steady {:.1f} m/s wind.\n'.format(wind_speed))
                f.write('! Time	Wind	Wind	Vert.	Horiz.	Vert.	LinV	Gust\n')
                f.write('!	Speed	Dir	Speed	Shear	Shear	Shear	Speed\n')
                f.write('  0.0\t{:.1f}\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n'.format(wind_speed))
                f.write('  0.1\t{:.1f}\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n'.format(wind_speed))
                f.write('999.9\t{:.1f}\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n'.format(wind_speed))
                
        # create FAST files
        jr.writeFASTFiles(turb_dir,TName,wind_fname,wind_speed,
                   BlPitch0=BlPitch0,RotSpeed0=RotSpeed0,
                   wind_dir=wind_dir,fileID=fileID,TMax=TMax)
        # run FAST
        print('Processing wind speed {:.1f}'.format(wind_speed))
        FASTfname = TName+'_'+fileID
        os.system('FAST.exe '+FASTfname+'.fst')
                      
        # load blade pitch angle, rotor speed
        FAST = jr.ReadFASTFile(FASTfname+'.out')
        BldPitch = FAST['Data'][:,FAST['Fields'].index('BldPitch2')]
        RotSpeed = FAST['Data'][:,FAST['Fields'].index('RotSpeed')]
                          
        # calculate and save last value. throw warning if not steady-state.
        BldPitchSS = np.mean(BldPitch[-10:])
        RotSpeedSS = np.mean(RotSpeed[-10:])
        if ((np.abs(BldPitch[0.83*BldPitch.size]-BldPitchSS)/BldPitchSS >= 0.10) \
                or (np.std(BldPitch[-50:])/BldPitchSS >= 0.05)):
            warnings.warn('Blade pitch angle is not at steady-state')
        elif ((np.abs(RotSpeed[0.83*RotSpeed.size]-RotSpeedSS)/RotSpeedSS >= 0.10) \
                or (np.std(RotSpeed[-50:])/RotSpeedSS >= 0.05)):
            warnings.warn('Rotor speed is not at steady-state')
        else:
            row = np.array([wind_speed,BldPitchSS,RotSpeedSS])
            LUT = np.vstack((LUT,row.reshape((1,row.size))))
                      
        # set initial conditions for next round, delete FAST files
        BlPitch0  = BldPitchSS*np.ones(3)
        RotSpeed0 = RotSpeedSS
        os.system('del '+FASTfname+'.fst')
        os.system('del '+FASTfname+'.out')
        os.system('del '+FASTfname+'_AD.ipt')
        
    else:
        BlPitch0  = LUT[np.where(LUT[:,0]==wind_speed)[0],1]*np.ones(3)
        RotSpeed0 = float(LUT[np.where(LUT[:,0]==wind_speed)[0],2])
        print('Wind speed {:.1f} present in LUT - skipping'.format(wind_speed))
    
print('Simulations completed.')    
   
# rearrange to increasing wind speed
LUT = LUT[LUT[:,0].argsort()]
   
# save LUT   
mdict = {}
mdict['LUT'] = LUT    
scio.savemat(fLUTpath,mdict)
print('Look-up table saved.')