"""
Script to do the following:
    1) iterate through NREL 20-Hz .mat files;
    2) load and clean time series;
    3) save screened/cleaned data in metadata table;
    4) save processed metadata.
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import scipy.io as scio

# get path to base directory
basedir = jr.getBasedir('NREL')

# get metadata fields
fields = jr.metadataFields('NREL')

# process metadata
metadata = jr.makeNRELmetadata(basedir)

# save dictionary
outdict = {}
outdict['fields']   = fields
outdict['metadata'] = metadata

# save as .mat file
fname = 'NREL_metadata_py.mat'
scio.savemat(fname,outdict)