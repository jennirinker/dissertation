#! /usr/bin/env python
""" Process statistics for specified wind direcotry
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import os
import numpy as np
import scipy.io as scio
from scipy.stats import skew
import JR_Library.main as jr

# =========== inputs =========== 

# path to turbine directory
TurbNames = ['WP0.75A08V00','WP1.5A08V03','WP3.0A02V02','WP5.0A04V00']
BaseTurbDir = '\\\\monsoon-data\\Public\\JRinker\\fast_simulations\\' + \
                'FastDir\\SmallRun'

# list of statistics to calculate for each channel 
calc_stats = ['max','min','mean','std','skew',
              '5%','25%','50%','75%','95%']

# =========== should not need to change =========== 

# define number of statistics to calculate
n_stats = len(calc_stats)

for TurbName in TurbNames:
    
    TurbDir = os.path.join(BaseTurbDir,TurbName)

    # get list of FAST output files in the directory
    orig_dir = os.getcwd()
    os.chdir(TurbDir)
    files = [f for f in os.listdir('.') if \
                ( os.path.isfile(f) and f.endswith('out') )]
    n_files = len(files)
    
    print('\nProcessing {:d} files from {:s}\n'.format( \
            n_files,TurbDir))
    
    # manually load first file to get list of variables
    fast_dict = jr.ReadFASTFile(files[0])
    fields    = fast_dict['Fields']
    n_fields  = len(fields)
    
    # loop through files and process statistics
    fnames     = []
    proc_stats = np.empty((n_files,n_fields*n_stats))
    for i_f in range(n_files):
        fname = files[i_f]
        fnames.append(fname)
    
        print('Processing file {:s} [{:d}/{:d}]...'.format( \
                fname,i_f,n_files))
    
        # load FAST file
        fast_dict   = jr.ReadFASTFile(fname)
        fast_fields = fast_dict['Fields']
        fast_data   = fast_dict['Data']
    
        # skip if list of fields does not match first file
        if set(fields).symmetric_difference(fast_fields):
            print('  *** skipping file {:s} (incorrect field list)'.format(fname))
            fnames[i_f]       = np.nan
            proc_stats[i_f,:] = np.nan
    
        else:
    
            # throw away first 30 seconds
            data = fast_data[np.where(fast_data[:,fields.index('Time')] \
                    >= 30.0),:][0]
    
            # sort time histories in increasing order
            data = np.sort(data,axis=0)
    
            # loop through statistics
            for i_stat in range(len(calc_stats)):
                if ( calc_stats[i_stat] == 'min' ):
                    proc_stats[i_f,i_stat*n_fields:(i_stat+1)*n_fields] = \
                        data[0,:]
                elif ( calc_stats[i_stat] == 'max' ):
                    proc_stats[i_f,i_stat*n_fields:(i_stat+1)*n_fields] = \
                        data[-1,:]
                elif ( calc_stats[i_stat] == 'mean' ):
                    proc_stats[i_f,i_stat*n_fields:(i_stat+1)*n_fields] = \
                        np.mean(data,axis=0)
                elif ( calc_stats[i_stat] == 'std' ):
                    proc_stats[i_f,i_stat*n_fields:(i_stat+1)*n_fields] = \
                        np.std(data,axis=0)
                elif ( calc_stats[i_stat] == 'skew' ):
                    proc_stats[i_f,i_stat*n_fields:(i_stat+1)*n_fields] = \
                        skew(data,axis=0)
                elif ( '%' in calc_stats[i_stat] ):
                    i_q = float(calc_stats[i_stat].rstrip('%'))/100.*data.shape[0]
                    proc_stats[i_f,i_stat*n_fields:(i_stat+1)*n_fields] = \
                        data[i_q,:]
                else:
                    err_str = 'Statistic {:s}'.format(calc_stats[i_stat]) + \
                        ' is not coded.'
                    ValueError(err_str)
               
    # create dictionary and save in script directory
    os.chdir(orig_dir)
    out_fname = TurbName + '_stats.mat'
    out_dict = {}
    out_dict['calc_stats'] = calc_stats
    out_dict['proc_stats'] = proc_stats
    out_dict['fields']     = fields
    out_dict['fnames']     = fnames
    scio.savemat(out_fname,out_dict)
    
    print('\nDictionary saved to {:s}\n'.format(out_fname))
