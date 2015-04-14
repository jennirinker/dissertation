""" Script to process data from a given dataset
"""

def processdata(dataset):
    """
    """
    import os
    import nrel_lib.nrel_io as nio
    import numpy as np

    # define metadata filename
    mdfname = 'metadata_' + dataset + '.txt'

    # if metadata file exists, ask permission to overwrite
    if os.path.isfile(mdfname):
        cont = input('Metadata file exists. Overwrite? [1/0]: ')
        if not cont:
            print '\nScript aborted.\n'
            return

    # get metadata fields for that dataset
    fields = nio.metadataFields(dataset)

    print '\nProcessing metadata.\n'

    # create metadata table
    metadata = nio.makemetadata(dataset)

    # save text file
    headerText = ','.join(fields)
    np.savetxt(mdfname,metadata,delimiter=',',header = headerText)

    print '\nMetadata saved.'
    print '\nScript complete.\n'
    
    return

# ==============================================================================

# specify which dataset to process
dataset = 'NREL'

# process the dataset and save the metadata
processdata(dataset)



