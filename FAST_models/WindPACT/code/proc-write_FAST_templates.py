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
DISCONInp   = 0

# define turbine name
turb_names = ['WP0.75A08V00','WP1.5A08V03','WP3.0A02V02','WP5.0A04V00']
turb_name  = turb_names[3]

# specify the directory to write the files to
turb_dir = os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
        'dissertation\\FAST_models\\FAST7',turb_name)
        
# specify turbine version if necessary
#turb_dir += '_newGBR'
#turb_dir += '_stifftwr'
        

# load turbine model
fTDictName = os.path.join(turb_dir,'parameters',turb_name+'_Dict.dat')
with open(fTDictName,'r') as f:
    TurbDict = json.load(f)

# interpolate blade/tower structural parameters and aerodynamic properties
#BldInterp, ADInterp = jr.InterpolateRotorParams(TurbDict)


if BldModesInp:
        
    # set filenames
    fname_temp = 'Template_Modes_Blade.inp'
    fname_out  = turb_name + '_BldModes.inp'
    
    sys.stdout.write('Writing Modes input file ' + \
                            '(blade) to \"{:s}\"...'.format(fname_out))
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(turb_dir,'modes',fname_out)
    
    # write blade file
    jr.writeBldModes(fpath_temp,fpath_out,TurbDict)
    
    sys.stdout.write('done.\n')
  
if TwrModesInp:
    
    # set filenames
    fname_temp = 'Template_Modes_Tower.inp'
    fname_out  = turb_name + '_TwrModes.inp'
    
    sys.stdout.write('Writing Modes input file' + \
                    ' (tower) to \"{:s}\"...'.format(fname_out))
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(turb_dir,'modes',fname_out)
    
    # write file
    jr.writeTwrModes(fpath_temp,fpath_out,TurbDict)
    
    sys.stdout.write('done.\n')
              
if BladeInp:
    
    # set filenames
    fname_temp = 'Template_Blade.dat'
    fname_out  = turb_name + '_Blade.dat'
    
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
    fname_out  = turb_name + '_AD.ipt'
    
    sys.stdout.write('Writing AeroDyn input file to \"{:s}\"...'.format(fname_out))
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(turb_dir,'templates',fname_out)
    
    # write blade file
    jr.WriteAeroDynTemplate(fpath_temp,fpath_out,TurbDict)
    
    sys.stdout.write('done.\n')
                
if TowerInp:
    
    # set filenames
    fname_temp = 'Template_Tower.dat'
    fname_out  = turb_name + '_Tower.dat'
    
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
    fname_out  = turb_name + '.fst'
    
    sys.stdout.write('Writing FAST input file to \"{:s}\"...'.format(fname_out))
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(turb_dir,'templates',fname_out)
    
    # write FAST file
    jr.WriteFASTTemplate(fpath_temp,fpath_out,TurbDict)
    
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
                
if DISCONInp:
    
    # set filenames
    fname_temp = 'Template_DISCON_nosat.f90'
    fname_out  = 'DISCON.f90'
    
    sys.stdout.write('Writing DISCON.f90 to \"{:s}\"...'.format(fname_out))
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(turb_dir,fname_out)
    
    # write Pitch.ipt
    jr.writeDISCON(fpath_temp,fpath_out,TurbDict)
                
    sys.stdout.write('done.\n')
    
                