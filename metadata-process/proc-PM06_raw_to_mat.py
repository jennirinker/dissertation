"""
Convert raw Plaine Morte data files to 10-minute .mat files
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import os, csv
import datetime as dt
import numpy as np
import scipy.io as scio

# define dataset-specific paramters
rec_mins = 10                       # number of minutes per record
header_line = 0                     # line number for file header
fields_line = 1                     # line number for data fields
units_line  = 2                     # line number for units
data_line   = 4                     # line number for start of data
N           = jr.datasetSpecs('PM06')[0]

def timestr2datetime(time_str):
    """ Convert string time to datetime structure """
    import datetime as dt
    
    if len(time_str) < 20:
        fstr = '%Y-%m-%d %H:%M:%S'
    else:
        fstr = '%Y-%m-%d %H:%M:%S.%f'
    time_dt = dt.datetime.strptime(time_str,fstr)
    
    return time_dt
    
def PM2numpy(row_str,time_0):
    """ Convert row of strings to numpy array """
    import numpy as np
    
    time_str = timestr2datetime(row_str[0])
    row_np = np.empty(len(row_str))
    for i in range(len(row_str)):
        if (i == 0):
            dt = (time_str - time_0).total_seconds()
            row_np[i] = dt
        elif ('NAN' in row_str[i]):
            row_np[i] = 'NAN'
        else:
            row_np[i] = float(row_str[i])
            
    return row_np

# define base directory for data
baseraw = 'H:\\data\\plaine-morte_raw\\' + \
    'CM 2006\\Data'
basepro = 'H:\\data\\plaine-morte'
# get list of folders
folders = [fold for fold in os.listdir(baseraw) if '-' in fold]

# loop through list of folders, processing each raw file
for i_fold in range(len(folders)):
    
    # get list of raw data files
    dirpath = os.path.join(baseraw,folders[i_fold])
    raws = [fname for fname in os.listdir(dirpath) if 'TS_WND' in fname]
    
    # loop through time series
    for i_file in range(len(raws)):
        
        # get file path
        fpath = os.path.join(dirpath,raws[i_file])
        print('Processing {}'.format(fpath))
        
        # process file
        with open(fpath,'r') as f:
            content = csv.reader(f,delimiter=',',quotechar='"')
            i_line = 0
            for row in content:
                
                # load header data if applicable
                if i_line == header_line:
                    header = row
                elif i_line == fields_line:
                    fields = row
                elif i_line == units_line:
                    units = row
                    
                # initialize data fields if it's first run into data
                elif i_line == data_line:
                    
                    # get first time stamp, define filename
                    time_str = row[0]
                    time_dt = timestr2datetime(time_str)
                    
                    
                    # define starting and ending time ranges
                    time_0  = time_dt
                    time_f  = time_dt - dt.timedelta(minutes= \
                                (time_dt.minute%rec_mins)-rec_mins,
                            seconds=time_dt.second,microseconds=time_dt.microsecond)
        
                    # initialize data parameters
                    i_t  = 0
                    data = np.empty((N,len(fields)))
                    
                elif i_line >= 4:
                    
                    # get the timestamp
                    time_str = row[0]
                    time_dt = timestr2datetime(time_str)
                    
                    # if the time is greater than the final measurement time, reset
                    if (time_dt >= time_f):
                        
                        # save data if enough time stamps are found in time period
                        if (i_t > N-5):
                            
                            # get/create directory name
                            dpath = os.path.join(basepro,
                                                 time_0.strftime('%Y\\%m\\%d'))
                            if (not os.path.isdir(dpath)): os.makedirs(dpath)
                            fname  = time_0.strftime('%m_%d_%Y_%H%M') + \
                                    '_' + os.path.splitext(raws[i_file])[0]
                            fpath  = os.path.join(dpath,fname)
                            values = data[:i_t,:]
                            
                            # create output dictionary
                            outdict           = {}
                            outdict['header'] = header
                            outdict['fields'] = fields
                            outdict['units']  = units
                            outdict['values'] = values
                            
                            # save output
                            scio.savemat(fpath,outdict)
                            
                            # print update statement
                            print('  '+time_0.strftime('%Y\\%m\\%d'))
                            
                        # reset values
                        i_t = 0
                        data = np.empty((N,len(fields)))
                        time_0 = time_f
                        time_f += dt.timedelta(minutes=rec_mins)
                                    
                    # save data in time range in numpy array
                    if ((time_dt >= time_0) and (time_dt < time_f)):
                        
                        # convert CSV list to numpy array
                        row_np = PM2numpy(row,time_0)
                        
                        # save numpy 1D array in 2D array
                        data[i_t,:] = row_np
                        
                        # update file index
                        i_t += 1
                        
                i_line += 1
             
                
