"""
Create a python dictionary describing a turbine from a series of text files
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import os
import numpy as np
import json



# define turbine name
TName,turb_id = 'WP0.75A08V00','WP0.75A08V00_stifftwr'
#TName = 'WP1.5A08V03'
#TName = 'WP3.0A02V02'
#TName = 'WP5.0A04V00'

# define directory locations: template files, turbine files
temp_dir = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation' + \
            '\\FAST_models\\WindPACT\\code\\templates'
turb_dir = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation' + \
            '\\FAST_models\\FAST7\\' + turb_id

# flag to indicate if blade/tower mode files are present
BModes = 1
TModes = 1

# ==================== Create and save turbine dictionary =====================
TurbDict = jr.createTurbineDictionary(TName,turb_dir,BModes=BModes,TModes=TModes)
fTDictName = os.path.join(turb_dir,'parameters\\'+TName + '_Dict.dat')
with open(fTDictName,'w') as f:
    json.dump(TurbDict,f)
    
print('\nProcessing complete.')
    
    
    