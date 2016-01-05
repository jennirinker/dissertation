"""
write FAST files for porting to peregrine
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\nwtc_python_tools'
if (libpath not in sys.path): sys.path.append(libpath)

import jr_fast
import os,json

# base turbine directory
BaseDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation' + \
            '\\FAST_models\FAST7'

RunName   = 'BigRun2'
TurbNames = ['WP0.75A08V00','WP1.5A08V03',
              'WP3.0A02V02','WP5.0A04V00']

TmplDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\' + \
            'nwtc_python_tools\\templates'

# directories and filepaths (shouldn't need to change for demo)
    

for TurbName in TurbNames:
    
    TurbDir = os.path.join(BaseDir,TurbName)
    AeroDir = os.path.join(TurbDir,'AeroData')
    DictPath = os.path.join(TurbDir,'parameters','{:s}_Dict.dat'.format(TurbName))
    
    # =============== should not need to change below this line ===================
    with open(DictPath,'r') as f_dict:
        TurbDict = json.load(f_dict)
    
    # write templates for files that depend on wind file (FAST and AeroDyn)
    jr_fast.WriteFAST7Template(TurbDict,TmplDir,TurbDir,TurbDir)
    jr_fast.WriteAeroDynTemplate(TurbDict,TmplDir,
                                 TurbDir,AeroDir,TurbDir)
    
    # write wind-independent files (blade files, tower files, and pitch file)
    jr_fast.WriteBladeFiles(TurbDict,TmplDir,TurbDir)
    jr_fast.WriteTowerFile(TurbDict,TmplDir,TurbDir)
    if (TurbDict['PCMode'] == 1):
        jr_fast.WritePitchCntrl(TurbDict,TmplDir,TurbDir)
