""" Messing with loading/processing metadata """

import sys
sys.path.append('..')

import JR_Library.data_analysis as da
import numpy as np
import matplotlib.pyplot as plt

# metadata filename
fname = 'metadata_NREL.txt'

# extract the fields and metadata from the text file if they aren't already
if ('metadata' not in vars()):
    
    # load data
    fields, metadata = da.loadmetadata(fname)

    # screen data
    cleandata = da.screenmetadata(fname,'NREL')
    
    # remove nans
    nanrows = np.unique(np.nonzero(np.isnan(cleandata))[0]) # rows with nans
    cleandata = np.delete(cleandata,nanrows,axis=0)         # delete nan rows


# variable columns
htCol = fields.index('Height')
rhouCol = fields.index('Concentration_u')

# let's look at 15 m
height = 15
data = cleandata[np.where(cleandata[:,htCol] == height)]

# extract concentration parameters
rho_u = data[:,rhouCol]

# plot a histogram
plt.figure()
plt.clf
