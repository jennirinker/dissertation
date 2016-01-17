"""
Load turbine model dictionary and create templates of FAST input files
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\nwtc_python_tools'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import json, os, sys, jr_fast

# which files to generate (1 = modes, 2 = fast/ad)
TmplType = 2

# define turbine name
turb_names = ['WP0.75A08V00','WP1.5A08V03','WP3.0A02V02','WP5.0A04V00']
turb_name  = turb_names[0]

# specify the directory to write the files to
TurbDir = os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
        'dissertation\\FAST_models\\FAST7',turb_name)
        
# specify turbine version if necessary
#turb_dir += '_newGBR'
#turb_dir += '_stifftwr'
        
TmplDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
            'public\\nwtc_python_tools\\templates'

# load turbine model
fTDictName = os.path.join(TurbDir,'parameters',turb_name+'_Dict.dat')
with open(fTDictName,'r') as f:
    TurbDict = json.load(f)

if TmplType == 1:
        
    # set filenames
    fname_temp = 'Template_Modes_Blade.inp'
    fname_out  = turb_name + '_BldModes.inp'
    
    sys.stdout.write('Writing Modes input file ' + \
                            '(blade) to \"{:s}\"...'.format(fname_out))
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(TurbDir,'modes',fname_out)
    
    # write blade file
    BldIdx = 1
    jr.writeBldModes(fpath_temp,fpath_out,TurbDict,BldIdx)
    
    sys.stdout.write('done.\n')
  
    # set filenames
    fname_temp = 'Template_Modes_Tower.inp'
    fname_out  = turb_name + '_TwrModes.inp'
    
    sys.stdout.write('Writing Modes input file' + \
                    ' (tower) to \"{:s}\"...'.format(fname_out))
    
    # set filepaths
    fpath_temp = os.path.join('templates',fname_temp)
    fpath_out  = os.path.join(TurbDir,'modes',fname_out)
    
    # write file
    jr.writeTwrModes(fpath_temp,fpath_out,TurbDict)
    
    sys.stdout.write('done.\n')
              
elif TmplType == 2:
    
    FastADTmplDir = os.path.join(TurbDir,'templates')
    AeroDir       = os.path.join(TurbDir,'AeroData')
    
    # write templates for files that depend on wind file (FAST and AeroDyn)
    jr_fast.WriteFAST7Template(TurbDict,TmplDir,TurbDir,FastADTmplDir)
    jr_fast.WriteAeroDynTemplate(TurbDict,TmplDir,
                                 TurbDir,AeroDir,FastADTmplDir)
    
    # write wind-independent files (blade files, tower files, and pitch file)
    jr_fast.WriteBladeFiles(TurbDict,TmplDir,TurbDir)
    jr_fast.WriteTowerFile(TurbDict,TmplDir,TurbDir)
    if (TurbDict['PCMode'] == 1):
        jr_fast.WritePitchCntrl(TurbDict,TmplDir,TurbDir)
    
                