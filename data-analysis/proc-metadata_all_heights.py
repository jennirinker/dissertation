"""
Load entire metadata set, screen it, isolate records with healthy data at all
heights, and save the result.
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import os
import numpy as np
import scipy.io as scio

# choose which dataset to process
#dataset = 'NREL'
#dataset = 'fluela'
#dataset = 'PM06'
dataset = 'texastech'

# define directory where wind parameters are stored (unused for matlab)
BaseDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                  'processed_data'

AllDataFields = ['Mean_Wind_Speed','Sigma_u','Tau_u','Concentration_u','Location_u',\
                  'Sigma_v','Tau_v','Concentration_v','Location_v',\
                  'Sigma_w','Tau_w','Concentration_w','Location_w','MO_Length_virt']
          
# ----------------------------------------------------------------------------

# import raw metadata
fields, clean = jr.loadmetadata(dataset)

# screen metadata, get measurement heights
screen   = jr.screenmetadata(fields,clean,dataset)
heights = jr.datasetSpecs(dataset)['IDs']

# loop through elements, searching to see if all heights are present
n_recs = screen.shape[0]
all_heights = np.empty((0,1+14*heights.size))
i_rec = 0
while i_rec < n_recs-heights.size:
    
    # if all height present
    if np.all(screen[i_rec:i_rec+heights.size,fields.index('ID')] == heights):
        
        # loop through heights, pulling out values and saving in row
        i_row = 0
        row = np.empty((1,1+len(AllDataFields)*heights.size))
        row[0,0] = screen[i_rec,0]
        for i_h in range(heights.size):
            
            for iparm in range(len(AllDataFields)):
                parm = AllDataFields[iparm]
                row[0,i_h*len(AllDataFields)+iparm+1] = \
                    screen[i_rec+i_h,fields.index(parm)]

        # append row to all rows
        all_heights = np.vstack((all_heights,row)) 
        
        # update loop index
        i_rec += heights.size
        
    # pass if not all heights present
    else:
        i_rec += 1
        
    # print update statement
    if not (i_rec % 100):
        print('{:.1f}'.format(i_rec/float(n_recs)*100))
             
        
NumOrig = len(screen)
NumAll  = len(all_heights)*len(heights)
print('Number of original records: {:d}'.format(NumOrig))
print('Number of original records: {:d} ({:.1f}%)'.format(NumAll,
                                      (NumOrig-NumAll)*100./NumOrig))
        
        
# save results in dictionary
outdict = {}
outdict['IDs']     = heights
outdict['all_fields']  = AllDataFields
outdict['values']  = all_heights
fname = dataset + '-metadata_allheights.mat'
fpath = os.path.join(BaseDir,fname)
scio.savemat(fpath,outdict)
print('\nData saved.')
