"""
Convert raw Plaine Morte data files to 10-minute .mat files
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import os, csv, shutil
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

# specify drive location
ext_drive = 'G:'

# flags to process intermediate, rewrite processed
proc_intr = 1
proc_proc = 1

# define base directory for intermediate data
baseraw = os.path.join(ext_drive,'data\\plaine-morte_raw\\' + \
    'CM 2006\\Data')
baseintr = os.path.join(ext_drive,'data\\plaine-morte_raw\\' + \
    'CM 2006\\intermediate_jr')
baseproc = os.path.join(ext_drive,'data\\plaine-morte\\PM06')

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

# process intermediate mat files
if proc_intr:
    
    # create intermediate directory if it doesn't exist
    if not os.path.isdir(baseintr):
        os.mkdir(baseintr)

    # get list of folders
    folders = [fold for fold in os.listdir(baseraw) if '-' in fold]
    
    # loop through list of folders, processing each raw file
    for i_fold in range(len(folders)):
        
        # get list of raw data files that don't end in 's'
        dirpath = os.path.join(baseraw,folders[i_fold])
        raws = [fname for fname in os.listdir(dirpath) if \
                (('TS_WND' in fname) and (not fname[:-4].endswith('s')))]
            
        # loop through files
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
                                dpath = os.path.join(baseintr,
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
   
# combine 3 separate mat files into single mat file         
if proc_proc:
    
    # delete processed directory if it exsits
    if os.path.isdir(baseproc):
        print('\nClearing processed directory...')
        shutil.rmtree(baseproc)
    os.mkdir(baseproc)
                
    # loop through year folders
    years = os.listdir(baseintr)
    for year in years:
        
        print('\nProcessing year {:s}'.format(year))
        
        # create year folder in processed directory
        year_dir = os.path.join(baseproc,year)
        os.mkdir(year_dir)
        
        # loop through months
        months = os.listdir(os.path.join(baseintr,year))
        for month in months:
            
            print('   month {:s}'.format(month))
            
            # create month folder in processed directory
            month_dir = os.path.join(year_dir,month)
            os.mkdir(month_dir)
            
            # loop through days
            days = os.listdir(os.path.join(baseintr,year,month))
            for day in days:
                
                print('      day {:s}'.format(day))
                
                # create day folder in processed directory
                day_dir = os.path.join(month_dir,day)
                os.mkdir(day_dir)
                
                # get path to intermediate files
                day_rdir = os.path.join(baseintr,year,month,day)
                
                # get list of unique file names without 1/2/3/2s/etc
                all_files = [fname[:-6] for fname in os.listdir(day_rdir) if \
                                (not fname[:-4].endswith('s'))]
                uniq_files = list(set(all_files))
                
                # loop through unique files
                for i_ufile in range(len(uniq_files)):
                    
                    uniq = uniq_files[i_ufile]
                    file_list = [fname for fname in os.listdir(day_rdir) if \
                                    uniq in fname]
                    
                    # try to combine all three mat files into one
                    # skip file on exception
                    try:
                        
                        # define path to save final dictionary
                        save_fname = uniq + '.mat'
                        save_fpath = os.path.join(day_dir,save_fname)
                        
                        # loop through all files
                        save_dict = {}
                        for i_file in range(len(file_list)):
                            ufname = file_list[i_file]
                            
                            # load intermediate dictionary
                            intr_dict = scio.loadmat(os.path.join(day_rdir,ufname))
                            intr_data = intr_dict['values']
                            intr_fields = [s.rstrip() for s in intr_dict['fields']]
                            
                            # save all fields individually in new dictionary
                            for i_field in range(len(intr_fields)):
                                field = intr_fields[i_field]
                                save_dict[field] = intr_data[:,i_field]
                        
                        # save dictionary
                        scio.savemat(save_fpath,save_dict)
                        
                    except IOError:
                        print('*** Error trying to load files for' + \
                                ' {:s}. Skipping.'.format(uniq))
                                
    print('\nProcessing complete\n')
            
            
        
    
    
