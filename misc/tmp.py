"""
testing routine to loop through M4 files, try to load them,
save a list of unloadable files, then re-download those files and try again
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import os
import scipy.io as scio
import numpy as np

# gat base directory for data
basedir = jr.getBasedir('NREL')

# define base URL
baseURL = 'http://wind.nrel.gov/MetData/' + \
    '135mData/M4Twr/20Hz/mat'
    
# walk through directory, try to load file, and save file path if can't
bad_files = []
i_files = 0
for root, dirs, files in os.walk(basedir,topdown=True):
    for f in files:
        fpath = os.path.join(root,f)
        if f.endswith('.mat'):
            try:
                struc = scio.loadmat(fpath)
                test = np.squeeze(struc['Sonic_u_15m'][0,0][0])
            except:
                print('Error loading {}'.format(fpath))
                bad_files.append(fpath)
            i_files += 1
            if (not i_files % 20):
                print(f)