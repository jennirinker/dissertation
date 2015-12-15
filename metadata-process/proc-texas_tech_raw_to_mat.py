"""
Process raw Texas Tech data to 10-minute .mat files
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import os, csv, re, shutil
import numpy as np
import scipy.io as scio

# define dataset-specific paramters
dataset     = 'texastech'
specs       =  jr.datasetSpecs(dataset)
n_t, dt, heights = specs['n_t'], specs['dt'], specs['IDs']

# specify drive location
ext_drive = 'G:'

# define base directory for intermediate data
baseraw  = ext_drive + '\\data\\texas-tech_raw'
baseproc = ext_drive + '\\data\\texas-tech'

# conversion factors
mph2mps     = 1./2.2369362920544    # mph to meters/s
vdc2mps_UVW = 40.26*mph2mps         # voltage to meters/s prop. anemometer
vdc2rh      = 20.                   # voltage to relative humidity
vdc2c_m     = 20.                   # voltage to celsius: slope
vdc2c_b     = -50.                  # voltage to celsius: offset
vdc2pscl_m  = 1.78*3386.39          # voltage to pascals: slope
vdc2pscl_b  = 23.6*3386.39          # voltage to pascals: offset

# create 30-minute time array
t_30min = np.arange(1800/dt)*dt

# delete processed directory if necessary
print('\nCleaning/creating previous processed directory...')
if os.path.isdir(baseproc):
    shutil.rmtree(baseproc)
os.makedirs(baseproc)

# get raw data field names
raw_field_fpath = 'G:\\data\\texas-tech_raw\\Documentation\\raw_fieldnames.csv'
with open(raw_field_fpath,'r') as csvfile:
    csvreader = csv.reader(csvfile)
    raw_fields_all = [[],[],[],[]]
    for row in csvreader:
        for i_row in range(4):
            raw_fields_all[i_row].append(row[i_row])
    for i_row in range(4):
        raw_fields_all[i_row] = [s for s in raw_fields_all[i_row] if s]
        
print('\nProcessing files...')

# for each set of 30-minute .csv files in raw data directory
data_folders = [fold for fold in os.listdir(baseraw) if fold[0].isdigit()]
for i_folder in range(len(data_folders)):
    data_folder = data_folders[i_folder]
    data_dir    = os.path.join(baseraw,data_folder)
        
    # get list of sets of 30-minute files
    files_all   = os.listdir(data_dir)
    files_trunc = [f[:-8] for f in files_all]
    files_uniq  = list(set(files_trunc))
    
    # loop through sets of 30-minite files
    for i_uniq in range(len(files_uniq)):
        
        uniq_fname = files_uniq[i_uniq]
        csv_list   = [fname for fname in os.listdir(data_dir) \
                        if uniq_fname in fname]

        # initialize 10-minute dictionaries (3 for 30-minute records) and
        # filenames for each mat file
        time_str0 = re.search('(_T[0-9]{4}_)',uniq_fname).group()[1:-1]
        time_str1 = time_str0[:3] + '1' + time_str0[-1]
        time_str2 = time_str0[:3] + '2' + time_str0[-1]
        fname0 = uniq_fname + '.mat'
        fname1 = re.sub('(T[0-9]{4})',time_str1,fname0)
        fname2 = re.sub('(T[0-9]{4})',time_str2,fname0)
        dict_10mins = [{'name':fname0},{'name':fname1},{'name':fname2}]
        
        # for raw files 1 through 4
        for i_f in range(len(csv_list)):
            
            # get filename and path to file
            csv_fname = csv_list[i_f]
            csv_fpath = os.path.join(data_dir,csv_fname)
            
            # extract list of raw fieldnames from list of all fieldnames
            raw_fields = raw_fields_all[i_f]
            int_fields = [s.replace('raw','int') for s in raw_fields]
            
            # process time series
            raw_data = np.genfromtxt(csv_fpath,delimiter=',')
            int_data = np.empty((t_30min.size,raw_data.shape[1]))
            nt_raw   = raw_data.shape[0]
            t_raw    = np.arange(nt_raw)/float(nt_raw)*1800
            for i_field in range(len(raw_fields)):
                field = raw_fields[i_field]
                
                # interpolate to even 30-minute base
                int_data[:,i_field] = np.interp(t_30min,
                                                t_raw,raw_data[:,i_field])
                                        
                # convert voltages to engineering units
                if ('UVW' in field):
                    int_data[:,i_field] = int_data[:,i_field] * vdc2mps_UVW
                elif ('Air_Temp' in field):
                    int_data[:,i_field] = int_data[:,i_field] * vdc2c_m + vdc2c_b
                elif ('Relative_Humidity' in field):
                    int_data[:,i_field] = int_data[:,i_field] * vdc2rh
                elif ('Baro_Presr' in field):
                    int_data[:,i_field] = int_data[:,i_field] * vdc2pscl_m + vdc2pscl_b
                    
            # split interpolated 30-min data into 10-minute sections and save
#            n_raw_10min = nt_raw/3.
            n_int_10min = 600./dt
            for i_t in range(3):
#                raw_data_10min = raw_data[(i_t*n_raw_10min):(i_t+1)*n_raw_10min]
                int_data_10min = int_data[(i_t*n_int_10min):(i_t+1)*n_int_10min]
                for i_field in range(len(raw_fields)):
                    field = raw_fields[i_field]
                    field_int = field.replace('_raw','')
                    if ('unused' not in field):
#                        dict_10mins[i_t][field]     = raw_data_10min[:,i_field]
                        dict_10mins[i_t][field_int] = int_data_10min[:,i_field]
                        
#                # rotate and save sonic anemometer data
#                sonic_x_keys = [key for key in int_fields if 'Sonic_x' in key]
#                for i_key in range(len(sonic_x_keys)):
#                    
#                    # get dictionary keys to sonic data
#                    sonic_x_key = sonic_x_keys[i_key]
#                    sonic_y_key = sonic_x_key.replace('Sonic_x','Sonic_y')
#                    sonic_z_key = sonic_x_key.replace('Sonic_x','Sonic_z')
#                    sonic_T_key = sonic_x_key.replace('Sonic_x','Sonic_T')
#                    
#                    # get raw sonic data
#                    x = dict_10mins[i_t][sonic_x_key]
#                    y = dict_10mins[i_t][sonic_y_key]
#                    z = dict_10mins[i_t][sonic_z_key]
#                    
#                    # rotate data
#                    u, v, w = jr.RotateTimeSeries(x,y,z)
#                    
#                    # save rotated data
#                    sonicrot_x_key = sonic_x_key.replace('x_int','u')
#                    sonicrot_y_key = sonic_y_key.replace('y_int','v')
#                    sonicrot_z_key = sonic_z_key.replace('z_int','w')
#                    sonicrot_T_key = sonic_T_key.replace('T_int','T')
#                    dict_10mins[i_t][sonicrot_x_key] = u
#                    dict_10mins[i_t][sonicrot_y_key] = v
#                    dict_10mins[i_t][sonicrot_z_key] = w
#                    dict_10mins[i_t][sonicrot_T_key] = dict_10mins[i_t][sonic_T_key]
#                    
#                # rotate and save propeller anemometer data
#                UVW_x_keys = [key for key in int_fields if 'UVW_x' in key]
#                for i_key in range(len(UVW_x_keys)):
#                    
#                    # get dictionary keys to propeller anemometer data
#                    UVW_x_key = UVW_x_keys[i_key]
#                    UVW_y_key = UVW_x_key.replace('UVW_x','UVW_y')
#                    UVW_z_key = UVW_x_key.replace('UVW_x','UVW_z')
#                    
#                    # get raw propeller anemometer data
#                    x = dict_10mins[i_t][UVW_x_key]
#                    y = dict_10mins[i_t][UVW_y_key]
#                    z = dict_10mins[i_t][UVW_z_key]
#                    
#                    # rotate data
#                    u, v, w = jr.RotateTimeSeries(x,y,z)
#                    
#                    # save rotated data
#                    UVWrot_x_key = UVW_x_key.replace('x_int','u')
#                    UVWrot_y_key = UVW_y_key.replace('y_int','v')
#                    UVWrot_z_key = UVW_z_key.replace('z_int','w')
#                    dict_10mins[i_t][UVWrot_x_key] = u
#                    dict_10mins[i_t][UVWrot_y_key] = v
#                    dict_10mins[i_t][UVWrot_z_key] = w
                    
        # get path to processed directory, create it if it doesn't exist
        date_str = re.search('(D[0-9]{8})',uniq_fname).group()[1:]
        year     = date_str[:4]
        month    = date_str[4:6]
        day      = date_str[6:8]
        proc_dir = os.path.join(baseproc,year,month,day)
        if not os.path.isdir(proc_dir):
            os.makedirs(proc_dir)
            print('  {:s}'.format(proc_dir))
                    
        # save 10-minute high-frequency dictionaries to .mat
        for i_t in range(3):
            fpath = os.path.join(proc_dir,dict_10mins[i_t]['name'])
            scio.savemat(fpath,dict_10mins[i_t])




    