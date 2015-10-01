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
siguCol = fields.index('Sigma_u')
tauuCol = fields.index('tau_u')
rhouCol = fields.index('Concentration_u')
muuCol  = fields.index('Location_u')
sigvCol = fields.index('Sigma_v')
tauvCol = fields.index('tau_v')
rhovCol = fields.index('Concentration_v')
muvCol  = fields.index('Location_v')
sigwCol = fields.index('Sigma_w')
tauwCol = fields.index('tau_w')
rhowCol = fields.index('Concentration_w')
muwCol  = fields.index('Location_w')
LMOCol = fields.index('MO_Length_virt')

# screen metadata, get measurement heights
clean   = jr.screenmetadata(fields,raw_parms,dataset)
heights = jr.datasetSpecs(dataset)[2]

# loop through elements, searching to see if all heights are present
n_recs = clean.shape[0]
all_heights = np.empty((0,1+14*heights.size))
i_rec = 0
fields = ['U','zeta','sig_u','L_u','rho_u','mu_u',\
                        'sig_v','L_v','rho_v','mu_v',\
                        'sig_w','L_w','rho_w','mu_w']
while i_rec < n_recs-heights.size:
    
    # if all height present
    if np.all(clean[i_rec:i_rec+heights.size,htCol] == heights):
        
        # loop through heights, pulling out values and saving in row
        i_row = 0
        row = np.empty((1,1+14*heights.size))
        row[0,0] = clean[i_rec,0]
        for i_h in range(heights.size):
            
            # extract values
            z   = clean[i_rec+i_h,htCol]
            U   = clean[i_rec+i_h,UCol]
            sigu = clean[i_rec+i_h,siguCol]
            tauu = clean[i_rec+i_h,tauuCol]
            rhou = clean[i_rec+i_h,rhouCol]
            muu  = clean[i_rec+i_h,muuCol]
            sigv = clean[i_rec+i_h,sigvCol]
            tauv = clean[i_rec+i_h,tauvCol]
            rhov = clean[i_rec+i_h,rhovCol]
            muv  = clean[i_rec+i_h,muvCol]
            sigw = clean[i_rec+i_h,sigwCol]
            tauw = clean[i_rec+i_h,tauwCol]
            rhow = clean[i_rec+i_h,rhowCol]
            muw  = clean[i_rec+i_h,muwCol]
            LMO = clean[i_rec+i_h,LMOCol]
            Lu   = tauu*U
            Lv   = tauv*U
            Lw   = tauw*U
            zMO = z/LMO
            
            # save in row
            row[0,i_h*len(fields) + 1]  = U
            row[0,i_h*len(fields) + 2]  = zMO
            row[0,i_h*len(fields) + 3]  = sigu
            row[0,i_h*len(fields) + 4]  = Lu
            row[0,i_h*len(fields) + 5]  = rhou
            row[0,i_h*len(fields) + 6]  = muu
            row[0,i_h*len(fields) + 7]  = sigv
            row[0,i_h*len(fields) + 8]  = Lv
            row[0,i_h*len(fields) + 9]  = rhov
            row[0,i_h*len(fields) + 10] = muv
            row[0,i_h*len(fields) + 11] = sigw
            row[0,i_h*len(fields) + 12] = Lw
            row[0,i_h*len(fields) + 13] = rhow
            row[0,i_h*len(fields) + 14] = muw
            
            
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
outdict['fields']  = fields
outdict['values']  = all_heights
fname = dataset + '-metadata_allheights.mat'
fpath = 'C:\\Users\\jrinker\\Dropbox\\research\\processed_data\\' + \
            fname
scio.savemat(fpath,outdict)
