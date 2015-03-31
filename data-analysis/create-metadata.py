""" Descriptive comment should go here
"""

import scipy.io as sio
import os
import numpy as np

# directory and filename path
# NEED DOUBLE SLASH ON WINDOWS
dname = 'G:\\data\\nrel-20Hz\\2013\\04\\09'
fname = '04_09_2013_02_40_00_049.mat'

# load mat file as Python dictionary
fpath = os.path.join(dname,fname)
struc = sio.loadmat(fpath)

# get time history
u = struc['Sonic_u_100m'][0,0][0]
dt = 1./20.
N = u.size
t = np.arange(0,N)*dt


# extract parameters to array

