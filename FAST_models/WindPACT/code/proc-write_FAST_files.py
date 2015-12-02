"""
taking turbine template files, specifying wind files and initial rotor speeds
and blade pitch angles
"""
import sys,os
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr

# specify a file ID and total simulation time
#fileID = '00000'
TMax = 630

# ******** set directory and turbine name ********
turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\FAST7\\WP0.75A08V00','WP0.75A08V00'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP0.75A08V00_newGBR','WP0.75A08V00'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP0.75A08V00_stiffblds','WP0.75A08V00'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP1.5A08V03','WP1.5A08V03'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP3.0A02V02','WP3.0A02V02'
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP5.0A04V00','WP5.0A04V00'
    
# set turbine version if necessary
ext = '_equil'
turb_dir += ext
    
# ******** set directory for wind ********
wind_dir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\wind_files'
#u0 = 15.0
#wind_fname = 'NoShr_'+'{:.1f}'.format(u0).zfill(4)+'.wnd'
#wind_fname = 'Step_04.0_05.0.wnd'
#wind_fname = 'Step_12.0_13.0.wnd'
#wind_fname  = '5m_TCno.bts'
#wind_fname  = 'Harmonic_11.7_3.0_0.2.wnd'
#wind_dir,fileID = turb_dir,'24134'
wind_dir = os.path.join(wind_dir,TName)
fileID = ['24134', '24142', '42331', '91242'][0]
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
