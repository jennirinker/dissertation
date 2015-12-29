"""
Create a python dictionary describing a turbine from a series of text files
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\nwtc_python_tools'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import jr_fast
import os
import json



# define turbine name
turb_names = ['WP0.75A08V00','WP1.5A08V03','WP3.0A02V02','WP5.0A04V00']
turb_name  = turb_names[0]

# specify the directory to write the files to
turb_dir = os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\' + \
        'dissertation\\FAST_models\\FAST7',turb_name)
        
# specify turbine version if necessary
#turb_dir += '_newGBR'
#turb_dir += '_stifftwr'

# define directory locations: template files, turbine files
temp_dir = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation' + \
            '\\FAST_models\\WindPACT\\code\\templates'
DefFastPath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation\\' + \
                'FAST_models\\WindPACT\\code\\templates\\WP_Template.fst'

# flag to indicate if blade/tower mode files are present
BModes = 1
TModes = 1

# ==================== Create and save turbine dictionary =====================

# load default dictionary
TurbDict = jr_fast.CreateFAST7Dict(DefFastPath)

# load custom parameters
NewDict = jr.CreateTurbineDictionary(turb_name,turb_dir,
                                      BModes=BModes,TModes=TModes)
                                      
# overwrite parameters in default dictionary
for key in NewDict.keys():
    if (key in TurbDict.keys()):
        TurbDict[key] = NewDict[key]
    else:
        ValueError('Wrong key {:s}'.format(key))
                                 
fTDictName = os.path.join(turb_dir,'parameters\\'+turb_name + '_Dict.dat')
with open(fTDictName,'w') as f:
    json.dump(TurbDict,f)
    
print('\nProcessing complete.')
    
    
    