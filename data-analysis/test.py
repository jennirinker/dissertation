""" Messing with loading/processing metadata """

import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
import JR_Library.main as jr

import numpy as np
import scipy.io as scio
# %%============================= load data ===================================

# path to matlab-processed metadata table
matpath = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
    '2015-02-28_temporal coherence in data\\code\\tempStruc.mat'

# load the structure
struc = scio.loadmat(matpath)

# extract the information
pdfTable = struc['tempStruc'][0,0][3]

iH = 0
iP = 0

# extract data
x = pdfTable[:,5*iH + iP]
x = np.sort(x)

fit_parms = jr.fitcompositedistribution('NREL',iP,x)
