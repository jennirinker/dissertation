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
raw_dir = 'G:\\data\\texas-tech_raw\\06-2012'

# time data
n_t, dt, heights = jr.datasetSpecs(dataset)
t_30min = np.arange(1800/dt)*dt

# raw data field names
raw_field_fpath = 'G:\\data\\texas-tech_raw\\Documentation\\raw_fieldnames.txt'
with open(raw_field_fpath,'r') as f:
    raw_fields_all = json.load(f)

# for each set of 30-minute .csv files in raw data directory
uniq_fname = 'FT2_E05_C01_R06405_D20120601_T1230_TR'
csv_list   = [fname for fname in os.listdir(raw_dir) if uniq_fname in fname]

# initialize 10-minute dictionaries (3 for 30-minute records)
dict_10mins = [{},{},{}]

# for raw files 1 through 4
for i_f in [3]:
    
    # get filename and path to file
    csv_fname = csv_list[i_f]
    csv_fpath = os.path.join(raw_dir,csv_fname)
    
    # extract list of raw fieldnames from list of all fieldnames
    raw_fields = raw_fields_all[i_f]
    int_fields = [s.replace('raw','int') for s in raw_fields]
    
    # interpolate all raw data to 30 min exactly (even unused data)
#    raw_data = np.genfromtxt(csv_fpath,delimiter=',')
    int_data = np.empty((t_30min.size,raw_data.shape[1]))
    nt_raw   = raw_data.shape[0]
    t_raw    = np.arange(nt_raw)/float(nt_raw)*1800
    for i_field in range(len(raw_fields)):
        field = raw_fields[i_field]
        int_data[:,i_field] = np.interp(t_30min,
                                        t_raw,raw_data[:,i_field])
    
    # split 30-min data into 10-minute sections
    n_raw_10min = nt_raw/3.
    n_int_10min = 600./dt
    
    # save used raw and interpolated data
    for i_t in range(3):
        raw_data_10min = raw_data[(i_t*n_raw_10min):(i_t+1)*n_raw_10min]
        int_data_10min = int_data[(i_t*n_int_10min):(i_t+1)*n_int_10min]
        for i_field in range(len(raw_fields)):
            field = raw_fields[i_field]
            field_int = field.replace('raw','int')
            if ('unused' not in field):
                dict_10mins[i_t][field]     = raw_data_10min[:,i_field]
                dict_10mins[i_t][field_int] = int_data_10min[:,i_field]
                
        # rotate and save sonic anemometer data
        sonic_x_keys = [key for key in int_fields if 'Sonic_x' in key]
        for i_key in range(len(sonic_x_keys)):
            
            # get dictionary keys to sonic data
            sonic_x_key = sonic_x_keys[i_key]
            sonic_y_key = sonic_x_key.replace('Sonic_x','Sonic_y')
            sonic_z_key = sonic_x_key.replace('Sonic_x','Sonic_z')
            
            # get raw sonic data
            x = dict_10mins[i_t][sonic_x_key]
            y = dict_10mins[i_t][sonic_y_key]
            z = dict_10mins[i_t][sonic_z_key]
            
            # rotate data
            x_raw = np.concatenate((x.reshape(x.size,1),
                          y.reshape(x.size,1),
                          z.reshape(x.size,1)),axis=1)
            x_rot, x_yaw = jr.RotateTimeSeries(x_raw)
            
            # save rotated data
            sonicrot_x_key = sonic_x_key.replace('x_int','u')
            sonicrot_y_key = sonic_y_key.replace('y_int','v')
            sonicrot_z_key = sonic_z_key.replace('z_int','w')
            dict_10mins[i_t][sonicrot_x_key] = x_rot[:,0]
            dict_10mins[i_t][sonicrot_y_key] = x_rot[:,1]
            dict_10mins[i_t][sonicrot_z_key] = x_rot[:,2]
            
        # rotate and save propeller anemometer data

    # convert voltages to engineering units

    