"""
Process raw Texas Tech data to 10-minute .mat files
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import os, json
import numpy as np


dataset = 'texastech'
raw_dir = 'G:\\data\\texas-tech_raw\\01-2012'

# time data
n_t, dt, heights = jr.datasetSpecs(dataset)
t_30min = np.arange(1800/dt)*dt

# raw data field names
raw_field_fpath = 'G:\\data\\texas-tech_raw\\Documentation\\raw_fieldnames.txt'
with open(raw_field_fpath,'r') as f:
    raw_fields_all = json.load(f)

# for each set of 30-minute .csv files in raw data directory
uniq_fname = 'FT2_E05_C01_R00001_D20120119_T2330_TR'
csv_list   = [fname for fname in os.listdir(raw_dir) if uniq_fname in fname]

# for raw files 1 through 4
for i_f in range(1):
    csv_fname = csv_list[i_f]
    csv_fpath = os.path.join(raw_dir,csv_fname)
    
    raw_fields = raw_fields_all[i_f]
    
    # interpolate raw data to 30 min exactly
#    raw_data = np.genfromtxt(csv_fpath,delimiter=',')
    int_data = np.empty((t_30min.size,raw_data.shape[1]))
    nt_raw   = raw_data.shape[0]
    t_raw    = np.arange(nt_raw)/float(nt_raw)*1800
    for i_field in range(len(raw_fields)):
        field = raw_fields[i_field]
        int_data[:,i_field] = np.interp(t_30min,
                                        t_raw,raw_data[:,i_field])
    
    # for each 10-minute section
    

    # rotate sonic data

    # convert voltages to engineering units

    