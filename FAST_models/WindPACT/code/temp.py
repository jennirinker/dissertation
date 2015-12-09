""" processing turbine models more automated
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import json, os

# define turbine name
turb_names = ['WP0.75A08V00','WP1.5A08V03','WP3.0A02V02''WP5.0A04V00']
turb_name  = turb_names[0]

# define directory locations: template files, turbine files
turb_dir = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation' + \
            '\\FAST_models\\FAST7\\' + turb_name

test_dict = 'Control'

# load saved turbine dictionary
fTDictName = os.path.join(turb_dir,'parameters',turb_name+'_Dict.dat')
with open(fTDictName,'r') as f:
    turb_dict_saved = json.load(f)

# compare with current turbine dictionary
turb_dict_new = jr.CreateTurbineDictionary(turb_name,turb_dir,
                                    BModes=0,TModes=0)

jr.WriteFASTTemplate(fpath_temp,fpath_out,TurbDict)


print('\nnew dictionary\n')

for key, value in sorted(turb_dict_new[test_dict].items()): # Note the () after items!
    print('{:20s} '.format(key),value)