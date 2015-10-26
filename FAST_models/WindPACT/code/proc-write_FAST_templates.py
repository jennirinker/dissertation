"""
Load turbine model dictionary and create templates of FAST input files
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import json, os, sys


BldModesInp = 0
TwrModesInp = 0
BladeInp    = 1
ADInp       = 1
TowerInp    = 1
FASTInp     = 1
PitchInp    = 1

# define turbine name
#TName = 'WP0.75A08V00'
#TName = 'WP1.5A08V03'
#TName = 'WP3.0A02V02'
TName = 'WP5.0A04V00'

# specify the directory to write the files to
turb_dir = os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
        'dissertation\\FAST_models\\FAST7',TName)

# load turbine model
fTDictName = os.path.join(turb_dir,'parameters',TName+'_Dict.dat')
with open(fTDictName,'r') as f:
    TurbDict = json.load(f)

# interpolate blade/tower structural parameters and aerodynamic properties
BldInterp, ADInterp = jr.InterpolateRotorParams(TurbDict)


if BldModesInp:
        
    # set filenames
    fname_temp = 'Template_Modes_Blade.inp'
    fname_out  = TName + '_BldModes.inp'
    
    sys.stdout.write('Writing Modes input file (blade) to \"{:s}\"...'.format(fname_out))
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(turb_dir,'modes',fname_out)
    
    # write blade file
    jr.writeBldModes(fpath_temp,fpath_out,TurbDict)
    
    sys.stdout.write('done.\n')
  
if TwrModesInp:
    
    # set filenames
    fname_temp = 'Template_Modes_Tower.inp'
    fname_out  = TName + '_TwrModes.inp'
    
    sys.stdout.write('Writing Modes input file (tower) to \"{:s}\"...'.format(fname_out))
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(turb_dir,'modes',fname_out)
    
    # write file
    jr.writeTwrModes(fpath_temp,fpath_out,TurbDict)
    
    sys.stdout.write('done.\n')
              
if BladeInp:
    
    # set filenames
    fname_temp = 'Template_Blade.dat'
    fname_out  = TName + '_Blade.dat'
    
    sys.stdout.write('Writing blade input file to \"{:s}\"...'.format(fname_out))
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(turb_dir,fname_out)
    
    # write blade file
    jr.writeBlade(fpath_temp,fpath_out,TurbDict)
    
    sys.stdout.write('done.\n')
                
if ADInp:
    
    # set filenames
    fname_temp = 'Template_AD.ipt'
    fname_out  = TName + '_AD.ipt'
    
    sys.stdout.write('Writing AeroDyn input file to \"{:s}\"...'.format(fname_out))
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(turb_dir,'templates',fname_out)
    
    # write blade file
    jr.writeAeroDynTemplate(fpath_temp,fpath_out,TurbDict)
    
    sys.stdout.write('done.\n')
                
if TowerInp:
    
    # set filenames
    fname_temp = 'Template_Tower.dat'
    fname_out  = TName + '_Tower.dat'
    
    sys.stdout.write('Writing tower input file to \"{:s}\"...'.format(fname_out))
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(turb_dir,fname_out)
    
    # write tower file
    jr.writeTower(fpath_temp,fpath_out,TurbDict)
    
    sys.stdout.write('done.\n')

if FASTInp:
    
    # set filenames
    fname_temp = 'Template.fst'
    fname_out  = TName + '.fst'
    
    sys.stdout.write('Writing FAST input file to \"{:s}\"...'.format(fname_out))
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(turb_dir,'templates',fname_out)
    
    # write FAST fi;e
    jr.writeFASTTemplate(fpath_temp,fpath_out,TurbDict)
    
    sys.stdout.write('done.\n')
                
if PitchInp:
    
    # set filenames
    fname_temp = 'Template_pitch.ipt'
    fname_out  = 'pitch.ipt'
    
    sys.stdout.write('Writing Pitch input file to \"{:s}\"...'.format(fname_out))
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(turb_dir,fname_out)
    
    # write Pitch.ipt
    jr.writePitch(fpath_temp,fpath_out,TurbDict)
                
    sys.stdout.write('done.\n')
                
                
                