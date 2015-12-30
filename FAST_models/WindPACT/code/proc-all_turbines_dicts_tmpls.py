"""
Re-process dictionaries and templates for all turbines
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\nwtc_python_tools'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import json, os, jr_fast


# define turbine name
TurbNames = ['WP0.75A08V00','WP1.5A08V03','WP3.0A02V02','WP5.0A04V00']

# specify the directory to write the files to
BaseDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
        'dissertation\\FAST_models\\FAST7'

DefFastPath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation\\' + \
                'FAST_models\\WindPACT\\code\\templates\\WP_Template.fst'
TmplDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
            'public\\nwtc_python_tools\\templates'
            

for TurbName in TurbNames:
    
    # ========================= create dictionary ===========================
    
    # load default dictionary
    TurbDict = jr_fast.CreateFAST7Dict(DefFastPath)
    
    # directory for that turbine
    TurbDir = os.path.join(BaseDir,TurbName)
    
    # load custom parameters
    NewDict = jr.CreateTurbineDictionary(TurbName,TurbDir)

    # overwrite parameters in default dictionary
    for key in NewDict.keys():
        TurbDict[key] = NewDict[key]
                                     
    fTDictName = os.path.join(TurbDir,'parameters\\'+ TurbName + '_Dict.dat')
    with open(fTDictName,'w') as f:
        json.dump(TurbDict,f)
        
    # ========================= write files ===========================
        
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
        
        