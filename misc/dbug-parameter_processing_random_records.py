"""
Testing the parallel wind parameter processing routine in Python
"""
import sys, os
from joblib import Parallel, delayed
import json
import numpy as np
import scipy.io as scio

# clean date 1: [2012,1,28,18,50]
# clean date 2: [2012,08,08,09,50]
# clean date 3: [2013,03,11,06,20]

if (__name__ == '__main__'):
    libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
    if (libpath not in sys.path): sys.path.append(libpath)
    import JR_Library.main as jr
    
    # variables for later
    basedir  = jr.getBasedir('NREL')
    njobs    = 4
    fields   = jr.metadataFields('NREL')
    n_recs   = 50
#    foutname = 'NREL-metadata.mat'
    
    # get indices of random mat files to proces
    lmats_fname = [fp for fp in os.listdir(basedir) if 'listmats' in fp][0]
    lmats_fpath = os.path.join(basedir,lmats_fname)
    with open(lmats_fpath,'r') as f:
        list_mats = json.load(f)
    n_files = len(list_mats)
    i_files = np.random.choice(n_files,n_recs,replace=False)
        
    # process files in parallel
    md_list = []
    for i in range(n_recs):
        md_list.append(jr.NRELlistmetadata(i_files[i],list_mats))
        print(i)
        sys.stdout.flush()
#    md_list = Parallel(n_jobs=njobs,verbose=9) \
#        (delayed(jr.NRELlistmetadata)(i_files[i],list_mats) for i in range(n_recs))
        
    # convert list of arrays to metadata
    metadata = np.empty((6*n_recs,len(fields)))
#    for i in range(n_files):
    for i in range(n_recs):
        metadata[6*i:6*(i+1)] = np.asarray(md_list[i])
        
#    # save output
#    mdict             = {}
#    mdict['fields']   = fields
#    mdict['metadata'] = metadata
#    scio.savemat(foutname,mdict)



# %% =============== draw samples of records, save metadata ===================
#n_recs  = 500                               # no. random records to draw
#yearRng = [2012,2015]                       # range of years
#fields = jr.metadataFields('NREL')          # MD fields
#
#CSLim  = 3                          # lower cup speed limit
#dir1   = 240                        # CCW edge for direction range
#dir2   = 315                        # CW edge for direction range
#preLim = 2.7                        # lower precipitation limit
#
# column indices for each value
#CScol  = fields.index('Wind_Speed_Cup')
#dirCol = fields.index('Wind_Direction')
#preCol = fields.index('Precipitation')
#
#i_rec = 0
#metadata = np.empty((n_recs,len(fields)))
#while i_rec < n_recs:
#    
#    # draw random time stamp/height
#    year   = random.randint(yearRng[0],yearRng[1])
#    month  = random.randint(1,12)
#    day    = random.randint(1,31)
#    hour   = random.randint(0,23)
#    minute = random.randint(0,6)*10
#    height = random.choice([15,30,50,76,100,131])
#    time_tup = (year,month,day,hour,minute)
#        
#    # check if it exists
#    fpath = jr.NRELtime2fpath(time_tup)
#    if len(fpath) > 0:
#        
#        # load structure
#        struc = scio.loadmat(fpath)
#        
#        # calculate parameters
#        row = jr.extractNRELparameters(struc,height)
#        
#        flagCS   = row[CScol] > CSLim
#        flagDir1 = row[dirCol] >= dir1
#        flagDir2 = row[dirCol] <= dir2
#        flagPrec = row[preCol] >= preLim
#        flagNaN  = not np.isnan(row[6])
#        
#        if (flagCS and flagDir1 and flagDir2 and flagPrec and flagNaN):
#            metadata[i_rec] = row.reshape(1,len(fields))
#    
#            i_rec += 1
#                                    
#            if not (i_rec % 5): print(i_rec)
                

