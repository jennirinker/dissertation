"""
Plot the wind velocity with the highest rho from the python-processed metadata
"""
import sys
sys.path.append('..')
import scipy.io as scio
import JR_Library.main as jr
import numpy as np

# load python-processed NREL metadata
fpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation\\' + \
    'nrel-data-process\\metadata_NREL.mat'
tmp = scio.loadmat(fpath)
fields = [string.rstrip(' ') for string in tmp['fields'].tolist()]
metadata = tmp['metadata']

# clean metadata
cleandata = jr.screenmetadata(fields,metadata,'NREL')

# find highest rho
rho_idx = fields.index('Concentration_u')
maxrho_idx = np.argmax(cleandata[:,rho_idx])

# load time history from that rho
timestamp = cleandata[maxrho_idx,0]
fpath = jr.NRELtime2fpath(timestamp)
struc = scio.loadmat(fpath)
t, u = jr.loadtimeseries('NREL',timestamp,15)
