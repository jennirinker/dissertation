"""
Load turbine model dictionary and create input files
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import json, os
import numpy as np

ModesInp = 0
BladeInp = 0
ADInp    = 0
TowerInp = 1

# define turbine name
#TName = 'WP0.75A08V00'
TName = 'WP1.5A08V03'
#TName = 'WP3.0A02V02'
#TName = 'WP5.0A04V00'

# specify the directory to write the files to
write_dir = os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
        'dissertation\\FAST_models\\FAST7',TName)

# load turbine model
fTDictName = TName + '_Dict.txt'
with open(fTDictName,'r') as f:
    TurbDict = json.load(f)

# interpolate blade/tower structural parameters and aerodynamic properties
BldInterp, ADInterp = jr.InterpolateRotorParams(TurbDict)
TowerInterp         = jr.InterpolateTowerParams(TurbDict)

if ModesInp:
    
    # set filenames
    fname_temp = 'Template_Modes_Blade.inp'
    fname_out  = TName + '_BldModes.inp'
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(write_dir,fname_out)
    
    # calculate modes-specific vales
    RotorRad = TurbDict['Rotor']['RotDiam']/2.
    SSAngVel = TurbDict['Nacelle']['RatedTipSpeed']/2/np.pi/RotorRad*60.
    HubRad   = TurbDict['Rotor']['HubDiam']/2
    
    # open template file and file to write to
    with open(fpath_temp,'r') as f_temp:
        with open(fpath_out,'w') as f_write:
            i_line = 0
            for line in f_temp:
                if i_line == 1:
                    f_write.write(line.format(SSAngVel))
                elif i_line == 3:
                    f_write.write(line.format(RotorRad))
                elif i_line == 4:
                    f_write.write(line.format(HubRad))
                elif i_line == 8:
                    f_write.write(line.format(len(BldInterp)))
                elif i_line == 12:
                    for i_BlNode in range(len(BldInterp)):
                        row = [BldInterp[i_BlNode,0],0.,
                               BldInterp[i_BlNode,3],
                               BldInterp[i_BlNode,4],BldInterp[i_BlNode,5]]
                        f_write.write(line.format(*row))
                else:
                    f_write.write(line)
                i_line += 1
                
if BladeInp:
    
    # set filenames
    fname_temp = 'Template_Blade.dat'
    fname_out  = TName + '_Blade.dat'
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(write_dir,fname_out)
    
    # calculate blade-specific vales
    title_str = 'FAST 7 blade file for turbine {:s}'.format(TName)
    Damping   = TurbDict['Rotor']['Damping']
    RotModes  = np.array(TurbDict['Rotor']['RotModes'])
    
    # open template file and file to write to
    with open(fpath_temp,'r') as f_temp:
        with open(fpath_out,'w') as f_write:
            i_line = 0
            for line in f_temp:
                if i_line == 2:
                    f_write.write(line.format(title_str))
                elif i_line == 4:
                    f_write.write(line.format(len(BldInterp)))
                elif i_line == 6:
                    f_write.write(line.format(Damping[0]*100.))
                elif i_line == 7:
                    f_write.write(line.format(Damping[1]*100.))
                elif i_line == 8:
                    f_write.write(line.format(Damping[2]*100.))
                elif i_line == 18:
                    for i_BlNode in range(len(BldInterp)):
                        f_write.write(line.format(*BldInterp[i_BlNode,:]))
                elif ((i_line >= 20) and (i_line <= 34)):
                    f_write.write(line.format(
                        RotModes.reshape(RotModes.size)[i_line-20]))
                else:
                    f_write.write(line)
                i_line += 1
                
if ADInp:
    
    # set filenames
    fname_temp = 'Template_AD.ipt'
    fname_out  = TName + '_AD.ipt'
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(write_dir,fname_out)
    
    # calculate blade-specific vales
    title_str = 'FAST 7 AeroDyn file for turbine {:s}'.format(TName)
    Airfoils  = TurbDict['Rotor']['Airfoils']
    
    # open template file and file to write to
    with open(fpath_temp,'r') as f_temp:
        with open(fpath_out,'w') as f_write:
            i_line = 0
            for line in f_temp:
                if i_line == 0:
                    f_write.write(line.format(title_str))
                elif i_line == 17:
                    f_write.write(line.format(len(Airfoils)))
                elif i_line == 18:
                    AFPath = 'AeroData/' + Airfoils[0] + '.dat'
                    f_write.write(line.format(AFPath))
                elif i_line == 19:
                    for i_AF in range(1,len(Airfoils)):
                        AFPath = 'AeroData/' + Airfoils[i_AF] + '.dat'
                        f_write.write(line.format(AFPath))
                elif i_line == 20:
                    f_write.write(line.format(len(ADInterp)))
                elif i_line == 22:
                    for i_AD in range(len(ADInterp)):
                        f_write.write(line.format(*ADInterp[i_AD,:]))
                else:
                    f_write.write(line)
                i_line += 1    
                
if TowerInp:
    
    # set filenames
    fname_temp = 'Template_Tower.ipt'
    fname_out  = TName + '_Tower.dat'
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(write_dir,fname_out)
    
    # calculate tower-specific vales
    title_str = 'FAST 7 tower file for turbine {:s}'.format(TName)
    Damping   = TurbDict['Rotor']['Damping']
    RotModes  = np.array(TurbDict['Rotor']['RotModes'])
    
    # open template file and file to write to
    with open(fpath_temp,'r') as f_temp:
        with open(fpath_out,'w') as f_write:
            i_line = 0
            for line in f_temp:
                if i_line == 2:
                    f_write.write(line.format(title_str))
                elif i_line == 4:
                    f_write.write(line.format(len(BldInterp)))
                elif i_line == 6:
                    f_write.write(line.format(Damping[0]*100.))
                elif i_line == 7:
                    f_write.write(line.format(Damping[1]*100.))
                elif i_line == 8:
                    f_write.write(line.format(Damping[2]*100.))
                elif i_line == 18:
                    for i_BlNode in range(len(BldInterp)):
                        f_write.write(line.format(*BldInterp[i_BlNode,:]))
                elif ((i_line >= 20) and (i_line <= 34)):
                    f_write.write(line.format(
                        RotModes.reshape(RotModes.size)[i_line-20]))
                else:
                    f_write.write(line)
                i_line += 1   
                
                
                
                