"""
write FAST files for porting to peregrine
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\nwtc_python_tools'
if (libpath not in sys.path): sys.path.append(libpath)
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import jr_fast
import os,json

# base turbine directory
BaseDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation' + \
            '\\FAST_models\FAST7'

RunName   = 'Fine'
RunParms = jr.RunName2WindParms(RunName)
TurbNames = RunParms['TurbNames']

TmplDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\' + \
            'nwtc_python_tools\\templates'
BaseModlDir = '\\\\monsoon-data\\Public\\JRinker\\fast_simulations\\ModlDir'
BaseWrDir   = '\\\\monsoon-data\\Public\\JRinker\\fast_simulations\\FastDir'
AeroDir     = '\\\\monsoon-data\\Public\\JRinker\\fast_simulations\\AeroDir'

# directories and filepaths (shouldn't need to change for demo)
    

for TurbName in TurbNames:
    
    TurbDir = os.path.join(BaseDir,TurbName)
    ModlDir = os.path.join(BaseModlDir,TurbName)
    FastDir = os.path.join(BaseWrDir,RunName,TurbName)
    DictPath = os.path.join(TurbDir,'parameters','{:s}_Dict.dat'.format(TurbName))
    FastADTmpDir = os.path.join(ModlDir,'templates')
    
    # =============== should not need to change below this line ===================
    with open(DictPath,'r') as f_dict:
        TurbDict = json.load(f_dict)
    
    # write templates for files that depend on wind file (FAST and AeroDyn)
    jr_fast.WriteFAST7Template(TurbDict,TmplDir,ModlDir,FastADTmpDir)
    jr_fast.WriteAeroDynTemplate(TurbDict,TmplDir,ModlDir,AeroDir,FastADTmpDir)

    
    # write wind-independent files (blade files, tower files, and pitch file)
    jr_fast.WriteBladeFiles(TurbDict,TmplDir,ModlDir)
    jr_fast.WriteTowerFile(TurbDict,TmplDir,ModlDir)
    if (TurbDict['PCMode'] == 1):
        jr_fast.WritePitchCntrl(TurbDict,TmplDir,TurbDir)
