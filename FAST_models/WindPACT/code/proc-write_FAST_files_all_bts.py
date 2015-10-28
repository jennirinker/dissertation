"""
Write FAST files for all .bts files in a directory
"""
import sys,os
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr

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

# specify total simulation time
TMax = 630
wind_dir = wr_dir

# get list of .bts files
bts_fnames = [f for f in os.listdir(wr_dir) if f.endswith('.bts')]
    
# write FAST files for each .bts file
for wind_fname in bts_fnames:
    jr.writeFASTFiles(turb_dir,turb_name,wind_fname,
                       wind_dir=wind_dir,TMax=TMax,
                       wr_dir=wr_dir)


