""" Script to process data from a given dataset
"""

def processdata(dataset):
    """
    """
    import os
    import nrel_lib.nrel_io as nio

    # define metadata filename
    mdfname = 'metadata_' + dataset + '.txt'

    # if metadata file exists, ask permission to overwrite
    if os.path.isfile(mdfname):
        cont = input('Metadata file exists. Overwrite? [1/0]: ')
        if not cont:
            return

    # get metadata fields for that dataset
    fields = nio.metadataFields(dataset)

    # create metadata table
    metadata = nio.makemetadata(dataset)

    print metadata

# ==============================================================================

# specify which dataset to process
dataset = 'NREL'

processdata(dataset)

# /media/jrinker/JRinker SeaGate External/data/nrel-20Hz/2012/02/13/02_13_2012_17_00_00_038.mat

