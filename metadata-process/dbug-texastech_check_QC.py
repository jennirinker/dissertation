"""
debugging texas tech metadata processing
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import scipy.io as scio
import os
import numpy as np

# define dataset
dataset = 'texastech'

BaseDataDir = 'C:\\Users\\jrinker\\Dropbox\\research\\processed_data'

# -----------------------------------------------------------------------------

# load metadata
DataName = '{:s}-metadata.mat'.format(dataset)
DataPath = os.path.join(BaseDataDir,DataName)
fields,metadata = jr.loadmetadata(DataPath)

# screen data by sonic speed and wind direction
cleandata = jr.screenmetadata(fields,metadata,dataset)

# extract largest wind speed/sig/etc
U = cleandata[:,fields.index('Mean_Wind_Speed')]
sig = cleandata[:,fields.index('Sigma_u')]
rho = cleandata[:,fields.index('Concentration_u')]

# get time series corresponding to max value(s)
idx = U.argmax()
ID = cleandata[idx,fields.index('ID')]
TSfield = 'Sonic_x'
RecTime = cleandata[idx,fields.index('Record_Time')]
TSdict = jr.loadtimeseries(dataset,TSfield,ID,RecTime)

