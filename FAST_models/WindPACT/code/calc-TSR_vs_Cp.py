"""
Plot TSR versus Cp

UNFINISHED*** couldn't figure out how to set rotor speed to constant
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
import os, json

# define TSRs, wind speed to analyze
TSRs = np.linspace(7,8,2)
wind_speed = 8.0
TMax = 10.0
rho  = 1.225

# set pitch angle to analyze
TName = 'WP1.5A08V03'

# specify the directory to write the files to
turb_dir = os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
        'dissertation\\FAST_models\\FAST7',TName)

# load turbine model
fTDictName = os.path.join(turb_dir,'parameters',TName+'_Dict.dat')
with open(fTDictName,'r') as f:
    TurbDict = json.load(f)
    
# set pitch angle to analyze
BldPitch = TurbDict['Rotor']['MinPitchAng']
RotorRad = TurbDict['Rotor']['RotDiam']/2.
GenRPM   = TurbDict['Nacelle']['RatedGenRPM']

# set wind speed directory, filename
wind_dir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\wind_files'
wind_fname = 'NoShr_'+'{:.1f}'.format(wind_speed).zfill(4)+'.wnd'
fileID     = ''

os.chdir(turb_dir)

# loop through TSRs
Cps = np.empty(TSRs.shape)
for i_u in range(len(TSRs)):
    
    # calculate rotor speed
    TSR      = TSRs[i_u]
    RotSpeed = TSR*wind_speed/RotorRad*(30./np.pi)
    
    # write FAST file
    jr.writeFASTFiles(turb_dir,TName,wind_fname,wind_speed,
                   wind_dir=wind_dir,fileID=fileID,TMax=TMax,GenDOF='False')
                   
    # run FAST
    fname = os.path.join(turb_dir,TName)
    os.system('FAST.exe '+fname+'.fst')
                   
    # load FAST file
    FAST = jr.ReadFASTFile(fname+'.out')
    
    # get generator power
    GenPwr = np.mean(FAST['Data'][:,FAST['Fields'].index('GenPwr')][-20:])
    
    # calculate and save cp value
    Cp = GenPwr/(0.5*rho*np.pi*(RotorRad**2)*wind_speed**3)
    Cps[i_u] = Cp
    
plt.figure(1)
plt.clf()
plt.plot(TSRs,Cps)