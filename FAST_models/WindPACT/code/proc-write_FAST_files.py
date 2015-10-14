"""
taking turbine template files, specifying wind files and initial rotor speeds
and blade pitch angles
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr

# specify a file ID and total simulation time
fileID = '00000'
TMax = 120

# set directories an filenames for wind
wind_dir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\wind_files'
u0 = 16.
wind_fname = 'NoShr_'+'{:.1f}'.format(u0).zfill(4)+'.wnd'
    
    
    
    
# ******** set directory and turbine name turbine 1 ********
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP1500_FAST_v7','WP1500'
#
## write files
#jr.writeFASTFiles(turb_dir,TName,wind_fname,u0,
#                   wind_dir=wind_dir,fileID=fileID,TMax=TMax)

# ******** set directory and turbine name turbine 2 ********
turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\FAST7\\WP1.5A08V03','WP1.5A08V03'
    
 # run FAST
jr.writeFASTFiles(turb_dir,TName,wind_fname,u0,
                   wind_dir=wind_dir,fileID=fileID,TMax=TMax)