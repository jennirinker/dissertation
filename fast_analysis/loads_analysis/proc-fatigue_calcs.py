"""
Debugging FAST fatigue calculations
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
import os
import numpy as np
import scipy.io as scio

# =========== inputs =========== 

# path to turbine directory
# path to turbine directory
TurbNames = ['WP0.75A08V00','WP1.5A08V03','WP3.0A02V02','WP5.0A04V00']
BaseTurbDir = '\\\\monsoon-data\\Public\\JRinker\\fast_simulations\\' + \
                'FastDir\\SmallRun'
#turb_dir = os.path.join('/scratch/jrinker/IOFiles',turb_name)
#FastDir = os.path.join('C:\\Users\\jrinker\\Documents\\GitHub\\dissertation\\' + \
#                'FAST_models\\FAST7',TurbName)

# if save dictionary
SaveDict = 1

# =========== should not need to change =========== 

for TurbName in TurbNames:
    
    FastDir = os.path.join(BaseTurbDir,TurbName)

    # get list of FAST output files in the directory
    files = [f for f in os.listdir(FastDir) if \
                (f.endswith('.out') )]
    n_files = len(files)
    
    print('\nProcessing {:d} files from {:s}\n'.format( \
            n_files,FastDir))
    
    # manually process first file to get list of DEL keys and slopes
    FastPath = os.path.join(FastDir,files[0])
    DELDict = jr.CalculateDELsAll(FastPath)
    DELs, DELkeys = DELDict['DELs'],DELDict['Keys']
    n_DELs  = len(DELkeys)
    
    # loop through files and process statistics
    fnames     = []
    DEL_stats = np.empty((n_files,n_DELs))
    for i_f in range(n_files):
        fname = files[i_f]
        fnames.append(fname)
    
        print('Processing file {:s} [{:d}/{:d}]...'.format( \
                fname,i_f+1,n_files))
    
        # process DEL info for FAST file
        FastPath = os.path.join(FastDir,fname)
        DELDict = jr.CalculateDELsAll(FastPath)
        fast_DELs, fast_DELkeys = DELDict['DELs'],DELDict['Keys']
    
        # skip if list of fields does not match first file
        if set(DELkeys).symmetric_difference(fast_DELkeys):
            print('  *** skipping file {:s} (incorrect field list)'.format(fname))
            fnames[i_f]       = np.nan
            DEL_stats[i_f,:] = np.nan
    
        else:
    
            # add data to DEL_stats
            for DEL_key in fast_DELkeys:
                DEL_stats[i_f,DELkeys.index(DEL_key)] = \
                    fast_DELs[fast_DELkeys.index(DEL_key)]
               
    # create dictionary and save in script directory
    if SaveDict:
        out_fname = TurbName + '_DELstats.mat'
        out_dict = {}
        out_dict['DELKeys'] = DELkeys
        out_dict['DELStats'] = DEL_stats
        out_dict['fnames']     = fnames
        scio.savemat(out_fname,out_dict)
        
        print('\nDictionary saved to {:s}\n'.format(out_fname))

