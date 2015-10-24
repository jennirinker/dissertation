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
#wind_dir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\wind_files'
u0 = 15.0
#wind_fname = 'NoShr_'+'{:.1f}'.format(u0).zfill(4)+'.wnd'
#wind_fname = 'Step_04.0_05.0.wnd'
#wind_fname = 'Step_12.0_13.0.wnd'
#wind_fname  = '5m_TCno.bts'
#wind_fname  = 'Kaimal_00000.wnd'


# ******** set directory and turbine name turbine 2 ********
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP0.75A08V00','WP0.75A08V00'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP1.5A08V03','WP1.5A08V03'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP3.0A02V02','WP3.0A02V02'
turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\FAST7\\WP5.0A04V00','WP5.0A04V00'
    
wind_dir = turb_dir
wind_fname = TName + '_' + fileID + '.bts'    
    
 # run FAST
jr.writeFASTFiles(turb_dir,TName,wind_fname,
                   wind_dir=wind_dir,fileID=fileID,TMax=TMax)

## ******** set directory and turbine name turbine 1 ********
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP1500_FAST_v7','WP1500'
#
## write files
#jr.writeFASTFiles(turb_dir,TName,wind_fname,
#                   wind_dir=wind_dir,fileID=fileID,TMax=TMax)
