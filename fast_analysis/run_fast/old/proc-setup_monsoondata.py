"""
Write Blade, Tower, pitch, and template files to monsoon
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\nwtc_python_tools'
if (libpath not in sys.path): sys.path.append(libpath)
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import jr_fast
import os, json


# define turbine name
TurbNames = ['WP0.75A08V00','WP1.5A08V03','WP3.0A02V02','WP5.0A04V00']

# directories and file paths
BaseDir = '\\\\monsoon-data\\Public\\JRinker\\fast_simulations\\ModlDir'
DefFastPath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation\\' + \
                'FAST_models\\WindPACT\\code\\templates\\WP_Template.fst'
TmplDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
            'public\\nwtc_python_tools\\templates'
AeroDir = '\\\\monsoon-data\\Public\\JRinker\\fast_simulations\\AeroDir'            
DictBaseDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation\\' + \
                'FAST_models\\FAST7'

# loop through 
for TurbName in TurbNames:
    
    # ========================= load dictionary ===========================
    
    # directory for that turbine
    TurbDir = os.path.join(BaseDir,TurbName)
    DictPath = os.path.join(DictBaseDir,TurbName,'parameters',
                            '{:s}_Dict.dat'.format(TurbName))
    with open(DictPath,'r') as fDict:
        TurbDict = json.load(fDict)
    
    # ========================= write files ===========================
        
    FastADTmplDir = os.path.join(TurbDir,'templates')
    
    # write templates for files that depend on wind file (FAST and AeroDyn)
    jr_fast.WriteFAST7Template(TurbDict,TmplDir,TurbDir,FastADTmplDir)
    jr_fast.WriteAeroDynTemplate(TurbDict,TmplDir,
                                 TurbDir,AeroDir,FastADTmplDir)
    
    # write wind-independent files (blade files, tower files, and pitch file)
    jr_fast.WriteBladeFiles(TurbDict,TmplDir,TurbDir)
    jr_fast.WriteTowerFile(TurbDict,TmplDir,TurbDir)
    if (TurbDict['PCMode'] == 1):
        jr_fast.WritePitchCntrl(TurbDict,TmplDir,TurbDir)
        

## define turbine name, turbine dictionary
## TurbName = ['WP0.75A08V00','WP1.5A08V03','WP3.0A02V02','WP5.0A04V00']
#TurbName = 'WP0.75A08V00'
#TurbDir = os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#        'dissertation\\FAST_models\\FAST7',TurbName)
#
## directories
#TmplDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\' + \
#            'nwtc_python_tools\\templates'     # loc of base templates
#ModlDir = os.path.join('\\\\monsoon-data\\Public\\JRinker\\' + \
#                    'fast_simulations\\ModlDir',TurbName)
#AeroDir = '\\\\monsoon-data\\Public\\JRinker\\' + \
#                    'fast_simulations\\AeroDir'
#TmplDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\' + \
#            'nwtc_python_tools\\templates'     # loc of base templates
##ModlDir = 'N:\\'
##AeroDir = 'R:\\'
##FastDir = 'A:\\'
##WindDir = 'P:\\'
##ModlDir_Wr       = os.path.join('\\\\monsoon-data\\Public\\JRinker\\' + \
##                    'fast_simulations\\ModlDir',TurbName)
#FastADTmplDir = os.path.join(ModlDir,'templates')
#WindDir = os.path.join('\\\\monsoon-data\\Public\\JRinker\\' + \
#                            'fast_simulations\\WindDir',TurbName)
#FastDir = os.path.join('\\\\monsoon-data\\Public\\JRinker\\' + \
#                        'fast_simulations\\FastDir',TurbName)
#
## -----------------------------------------------------------------------------
#
## load turbine dictionary from FAST file
#FastName = [f for f in os.listdir(TurbDir) if f.endswith('.fst')][0]
#TurbDict = jr_fast.CreateFAST7Dict(os.path.join(TurbDir,FastName))
#TurbDict['TurbName'] = TurbName
#
## write templates for files that depend on wind file (FAST and AeroDyn)
#jr_fast.WriteFAST7Template(TurbDict,TmplDir,ModlDir,FastADTmplDir)
#jr_fast.WriteAeroDynTemplate(TurbDict,TmplDir,
#                             ModlDir,AeroDir,FastADTmplDir)
#
## write wind-independent files (blade files, tower files, and pitch file)
#jr_fast.WriteBladeFiles(TurbDict,TmplDir,ModlDir)
#jr_fast.WriteTowerFile(TurbDict,TmplDir,ModlDir)
#if (TurbDict['PCMode'] == 1):
#    jr_fast.WritePitchCntrl(TurbDict,TmplDir,FastDir)
#    
#
## write wind-dependent files (FAST and AeroDyn files) for all wind files in
##   specified directory
#SimSpecs = {'TMax':630.,'TStart':30.}
#jr_fast.WriteFastADAll(TurbName,ModlDir,WindDir,FastDir,
#                            Naming=1,**SimSpecs)

