""" write wind files
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr

# set wind parameters
n_osc, A, U, f  = 3, 3, 11.7, 0.2
T, T_steady, dt = 630., 30., 0.05

# set file ID
fileID = '{:.1f}_{:.1f}_{:.1f}'.format(U,A,f)

# set wind directory
wind_dir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
                'dissertation\\FAST_models\\wind_files'

# write files
jr.writeHarmonicWind(U,A,f,T=T,dt=dt,T_steady=T_steady,
                    n_osc=4,fileID=fileID,wind_dir=wind_dir)