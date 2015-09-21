"""
Load entire metadata set, screen it, isolate records with clean data at all
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
dataset = 'NREL'

# define directory where wind parameters are stored (unused for matlab)
basedir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                  'processed_data'

# import raw metadata
fname = dataset + '-metadata.mat'
fpath = os.path.join(basedir,fname)
fields, raw_parms = jr.loadmetadata(fpath)

# get columns
htCol  = fields.index('Height')
UCol   = fields.index('Mean_Wind_Speed')
sigCol = fields.index('Sigma_u')
tauCol = fields.index('tau_u')
rhoCol = fields.index('Concentration_u')
muCol  = fields.index('Location_u')
LMOCol = fields.index('MO_Length_interp')

# screen metadata, get measurement heights
clean   = jr.screenmetadata(fields,raw_parms,dataset)
heights = jr.datasetSpecs(dataset)[2]

# loop through elements, searching to see if all heights are present
n_recs = clean.shape[0]
all_heights = np.empty((0,1+6*heights.size))
i_rec = 0
while i_rec < n_recs-heights.size:
    
    # if all height present
    if np.all(clean[i_rec:i_rec+heights.size,htCol] == heights):
        
        # loop through heights, pulling out values and saving in row
        i_row = 0
        row = np.empty((1,1+6*heights.size))
        row[0,0] = clean[i_rec,0]
        for i_h in range(heights.size):
            
            # extract values
            z   = clean[i_rec+i_h,htCol]
            U   = clean[i_rec+i_h,UCol]
            sig = clean[i_rec+i_h,sigCol]
            tau = clean[i_rec+i_h,tauCol]
            rho = clean[i_rec+i_h,rhoCol]
            mu  = clean[i_rec+i_h,muCol]
            LMO = clean[i_rec+i_h,LMOCol]
            L   = tau*U
            zMO = z/LMO
            
            # save in row
            row[0,i_h*heights.size + 1] = U
            row[0,i_h*heights.size + 2] = sig
            row[0,i_h*heights.size + 3] = L
            row[0,i_h*heights.size + 4] = rho
            row[0,i_h*heights.size + 5] = mu
            row[0,i_h*heights.size + 6] = zMO
            
        # append row to all rows
        all_heights = np.vstack((all_heights,row)) 
        
        # update loop index
        i_rec += heights.size
        
    # pass if not all heights present
    else:
        i_rec += 1
        
    # print update statement
    if not (i_rec % 100):
        print(i_rec/float(n_recs)*100)
        
# save results in dictionary
outdict = {}
outdict['heights'] = heights
outdict['fields']  = ['U','sig','L','rho','mu','zeta']
outdict['values']  = all_heights
fname = dataset + '-metadata_allheights.mat'
scio.savemat(fname,outdict)
