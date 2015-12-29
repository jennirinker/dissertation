"""
Write FAST files for all .bts files in a directory
"""
import sys,os
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\nwtc_python_tools'
if (libpath not in sys.path): sys.path.append(libpath)

#import JR_Library.main as jr
import jr_fast

# turbine names
turb_names = ['WP0.75A08V00','WP1.5A08V03','WP3.0A02V02','WP5.0A04V00']
TurbName = turb_names[1]
#TurbName = 'WP1500'

# specity turbine name, read directory, and write directory
base_wdir = 'C:\\Users\\jrinker\\Documents\\GitHub' + \
        '\\dissertation\\FAST_models\\verification'
wr_dir = os.path.join(base_wdir,TurbName)
base_rdir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\FAST7'
TurbDir = os.path.join(base_rdir,TurbName)

# specify turbine version if necessary
#turb_dir += '_newGBR'
#turb_dir += '_FAST_v7'

# specify wind directory
WindDir = os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\dissertation\\' + \
                'FAST_models\\wind_files',TurbName)
#wind_dir = os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\dissertation\\' + \
#                'FAST_models\\wind_files','WP1.5A08V03')
                
jr_fast.WriteFastADAll(TurbName,TurbDir,WindDir,TurbDir)




