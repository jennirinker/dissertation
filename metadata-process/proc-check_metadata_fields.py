"""
A script to run <calculatefield> and get the output dictionary fields to put
into <metadataFields>
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import scipy.io as scio

# pick dataset
#datasets = ['NREL','fluela','PM06','texastech']
dataset = 'PM06'

# load high-frequency structure
if dataset == 'NREL':
    fpath_hf = 'G:\\data\\nrel-20Hz\\2012\\03\\02\\03_02_2012_01_00_00_032.mat'
elif dataset == 'fluela':
    fpath_hf = 'G:\\data\\fluela-high_freq\\2010\\' + \
            '02\\05\\02_05_2010_13_30_00_000.mat'
elif dataset == 'PM06':
    fpath_hf = 'G:\\data\\plaine-morte\\CM06\\2006\\02\\' + \
                    '02\\02_02_2006_1440_TS_WND.mat'
elif dataset == 'texastech':
    fpath_hf = ''
struc_hf = scio.loadmat(fpath_hf)

print('\nComparing metadata fields for dataset {:s}'.format(dataset))

# calculate fields in output dictionary
IDs = jr.datasetSpecs(dataset)['IDs']
outdict = jr.calculatefield(dataset,struc_hf,IDs[-1])
dict_keys = outdict.keys()

# load fields stored in module
mod_keys = jr.metadataFields(dataset)

# load fields stored in metadata
fpath = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
            'processed_data\\{:s}-metadata.mat'.format(dataset)
fields, raw_md = jr.loadmetadata(fpath)
md_keys = fields

# print intersection
print('\nDifference between calculatfield and metadataFields:')
print('----------------------------------------------------')
print(set(dict_keys) ^ set(mod_keys))

print('\nDifference between stored fields and metadataFields:')
print('----------------------------------------------------')
print(set(md_keys) ^ set(mod_keys))
