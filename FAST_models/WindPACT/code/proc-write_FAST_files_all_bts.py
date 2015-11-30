"""
Write FAST files for all .bts files in a directory
"""
import sys,os
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

#import JR_Library.main as jr
import JR_Library.jr_fast as jr_fast

# turbine names
turb_names = ['WP0.75A08V00','WP1.5A08V03','WP3.0A02V02','WP5.0A04V00']
turb_name = turb_names[0]

# specity turbine name, read directory, and write directory
base_wdir = 'C:\\Users\\jrinker\\Documents\\GitHub' + \
        '\\dissertation\\FAST_models\\verification'
wr_dir = os.path.join(base_wdir,turb_name)
base_rdir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\FAST7'
turb_dir = os.path.join(base_rdir,turb_name)

# specify turbine version if necessary
#turb_dir += '_newGBR'

# specify wind directory
wind_dir = os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\dissertation\\' + \
                'FAST_models\\wind_files',turb_name)

jr_fast.WriteFAST7InputsAll(turb_dir,turb_name,wind_dir)



