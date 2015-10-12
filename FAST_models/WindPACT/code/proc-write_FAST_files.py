"""
taking turbine template files, specifying wind files and initial rotor speeds
and blade pitch angles
"""
import numpy as np
import os

# specify a file ID and total simulation time
fileID = '_00000'
TMax = 120

# set directories an filenames for wind
wind_dir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\wind_files'
wind_fname = 'NoShr_7.wnd'    
    
# set directory and turbine name
#turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP1500_FAST_v7','WP1500'
turb_dir,TName = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\FAST7\\WP1.5A08V03','WP1.5A08V03'

# set initial conditions
BlPitch0  = 2.6*np.ones(3)
RotSpeed0 = 20.0

# create filenames
fAD_name  = TName+fileID+'_AD.ipt'
fFST_name = TName+fileID+'.fst'
fAD_temp = os.path.join(turb_dir,'templates',TName+'_AD.ipt')
fAD_out    = os.path.join(turb_dir,fAD_name)
fFST_temp = os.path.join(turb_dir,'templates',TName+'.fst')
fFST_out  = os.path.join(turb_dir,fFST_name)

# write AeroDyn file
with open(fAD_temp,'r') as f_temp:
    with open(fAD_out,'w') as f_write:
        i_line = 0
        for line in f_temp:
            if i_line == 9:
                f_write.write(line.format(os.path.join(wind_dir,wind_fname)))
            else:
                f_write.write(line)
            i_line += 1
            
# write FAST file
with open(fFST_temp,'r') as f_temp:
    with open(fFST_out,'w') as f_write:
        i_line = 0
        for line in f_temp:
            if i_line == 9:
                f_write.write(line.format(TMax))
            elif i_line == 45:
                f_write.write(line.format(BlPitch0[0]))
            elif i_line == 46:
                f_write.write(line.format(BlPitch0[1]))
            elif i_line == 47:
                f_write.write(line.format(BlPitch0[2]))
            elif i_line == 72:
                f_write.write(line.format(RotSpeed0))
            elif i_line == 160:
                f_write.write(line.format(fAD_name))
            else:
                f_write.write(line)
            i_line += 1