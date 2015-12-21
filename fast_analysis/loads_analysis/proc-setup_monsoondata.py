"""
Set up model, wind, and aerodata directories in //monsoon-data
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\nwtc_python_tools'
if (libpath not in sys.path): sys.path.append(libpath)
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import jr_fast
import JR_Library.main as jr
import os

# define turbine name, turbine dictionary
# TurbName = ['WP0.75A08V00','WP1.5A08V03','WP3.0A02V02','WP5.0A04V00']
TurbName = 'WP0.75A08V00'
TurbDir = os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
        'dissertation\\FAST_models\\FAST7',TurbName)

# directories
TmplDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\' + \
            'nwtc_python_tools\\templates'     # loc of base templates
ModlDir = os.path.join('\\\\monsoon-data\\Public\\JRinker\\' + \
                    'fast_simulations\\ModlDir',TurbName)
AeroDir = '\\\\monsoon-data\\Public\\JRinker\\' + \
                    'fast_simulations\\AeroDir'
TmplDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\' + \
            'nwtc_python_tools\\templates'     # loc of base templates
#ModlDir = 'N:\\'
#AeroDir = 'R:\\'
#FastDir = 'A:\\'
#WindDir = 'P:\\'
#ModlDir_Wr       = os.path.join('\\\\monsoon-data\\Public\\JRinker\\' + \
#                    'fast_simulations\\ModlDir',TurbName)
#WrDir_FastADTmpl = os.path.join(ModlDir_Wr,'templates')
WindDir = os.path.join('\\\\monsoon-data\\Public\\JRinker\\' + \
                            'fast_simulations\\WindDir',TurbName)
FastDir = os.path.join('\\\\monsoon-data\\Public\\JRinker\\' + \
                        'fast_simulations\\FastDir',TurbName)

# -----------------------------------------------------------------------------

# load turbine dictionary from FAST file
FastName = [f for f in os.listdir(TurbDir) if f.endswith('.fst')][0]
TurbDict = jr_fast.CreateFAST7Dict(os.path.join(TurbDir,FastName))
TurbDict['TurbName'] = TurbName

# write templates for files that depend on wind file (FAST and AeroDyn)
jr_fast.WriteFAST7Template(TurbDict,TmplDir,ModlDir)
jr_fast.WriteAeroDynTemplate(TurbDict,TmplDir,
                             ModlDir,AeroDir)

# write wind-independent files (blade files, tower files, and pitch file)
jr_fast.WriteBladeFiles(TurbDict,TmplDir,ModlDir)
jr_fast.WriteTowerFile(TurbDict,TmplDir,ModlDir)
if (TurbDict['PCMode'] == 1):
    jr_fast.WritePitchCntrl(TurbDict,TmplDir,ModlDir)
    

# write wind-dependent files (FAST and AeroDyn files) for all wind files in
#   specified directory
SimSpecs = {'TMax':630.,'TStart':30.}
jr_fast.WriteFastADAll(TurbName,ModlDir,WindDir,FastDir,
                            Naming=1,**SimSpecs)

