"""
Library of functions used in dissertation analyses.

Module is separated into groups of functions:
    - File I/O:
        Loading metadata files, TurbSim files, time series
    - Metadata processing:
        Extracting wind parameters from data
    - Metadata analysis:
        Calculating NSAE, fitting composite distributions
    - Simulation:
        Kaimal simulation, IEC functions
    - TurbSim Analysis:
        Hub-height pase differences, PSDs, velocity profiles
    - Mappings:
        Time UTC to local, general field to dataset-specific, etc.
    - Miscellaneous:
        Other functions

Restructured 2015-06-30
Jenni Rinker, Duke University
"""

# %%===========================================================================
# FILE I/O
# =============================================================================

def loadmetadata(fname):
    """ Load data and fields from metadata file processed in Python,
        either .txt or .mat
    
        Args:
            fname (string): path to file
            
        Returns:
            fields (list): fields in each column
            metadata (numpy array): values for each field and each record
    """
    import numpy as np
    import scipy.io as scio

    # if it is a text file
    if (fname.endswith('txt')):
    
        # strip fields from header line
        with open(fname,'r') as f:
            fields = f.readline().lstrip('# ').rstrip('\n').split(',')
            
        # read in data table
        metadata = np.loadtxt(fname,delimiter=',',skiprows=1)

    # if it is a .mat file
    elif (fname.endswith('mat')):

        struc     = scio.loadmat(fname)         # load structure
        fields    = [string.rstrip(' ') for \
                  string in  struc['fields']]   # load fields, stripping spaces
        metadata  = struc['metadata']

    else:
        errStr = 'Invalid file type \"{}\"'.format(fname)
        raise TypeError(errStr)
    
    return (fields, metadata)


def loadNRELmatlab():
    """ Load the (uncleaned) metadata table processed in Matlab

        Returns:
            fields (list): names of fields in coumns
            metadata (numpy array): array of atmospheric params
    """
    import scipy.io as scio
    import numpy as np
    
    # path to matlab-processed metadata table
    matpath = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                  'processed_data\\NREL-mat-metadata.mat'

    # load the structure
    struc = scio.loadmat(matpath)

    # extract the information
    dataStruc = struc['dataStruc'][0,0]

    # metadata table
    metadata = dataStruc[0]

    # oddly-formatted field names
    tmp = np.squeeze(dataStruc[1])
        
    # extract readable fieldnames
    fields = [str(ele[0]) for ele in tmp]

    # replace necessary old fieldnames with new ones for screening
    fields[fields.index('Cup_speed')]      = 'Wind_Speed_Cup'
    fields[fields.index('Wind_direction')] = 'Wind_Direction'
    fields[fields.index('Precip_intens')]  = 'Precipitation'
    fields[fields.index('rho_u')]          = 'Concentration_u'
    fields[fields.index('U')]              = 'Mean_Wind_Speed'
    fields[fields.index('sigma_u')]        = 'Sigma_u'
    fields[fields.index('mu_u')]           = 'Location_u'

    return (fields, metadata)

def loadtimeseries(dataset,field,ht,data_in):
    """ Load the time series for a given dataset, field and timestamp 
        or structure

        Args:
            dataset (string): flag to indicate dataset
            field (string): general field name
            ht (int): measurement height
            data_in (float/tuple or dictionary): float or tuple representing 
                time value or high-frequency data structure

        Returns:
            outdict (dictionary): keys = ['raw','clean','flags']
    """
    import scipy.io as scio
    import numpy as np
    
    # make sure datfield is valid before proceeding
    datfield = field2datfield(dataset,field,ht)
    if (not check_datfields(dataset,datfield)):
        raise KeyError('Invalid datafield \"{}\".'.format(\
                    datfield))

    # turn off warnings about invalid entries
    np.seterr(invalid='ignore')
    
    # define max number of spikes
    n_spikes_max = 20

    if (dataset == 'NREL'):

        # set data ranges, desired length of time series
        dataRng = dataRanges(dataset,datfield)
        dt, N   = 0.05, 12000
        t       = np.arange(N)*dt
                    
        # determine what kind of data was provided
        if (isinstance(data_in,dict)):                   # data structure given
            struc = data_in
        elif (isinstance(data_in,(int,float,tuple))):    # timestamp given
            fpath = time2fpath(dataset,data_in)
	    # try to load the structure
	    try:
		struc = scio.loadmat(fpath)
	    except:
		errStr = 'Corrupt or empty structure for {}'.format(fpath)
		raise TypeError(errStr)
        else:
            errStr = 'Incorrect data type ' + \
                '\"{}\" for input \"data_in\".'.format(type(data_in))
            raise TypeError(errStr)

        # initialize list of flags
        flags = []

        # try to load time series
        try:
            x_raw = np.squeeze(struc[datfield][0,0][0])            
            # flag 5003: data are all NaNs
            if (np.all(np.isnan(x_raw))):
                flags.append(5003)                
        except:
            # flag 1006: data are not in 20-Hz structure
            x_raw = np.empty(N)
            x_raw[:] = np.nan
            flags.append(1006)
                
        # copy raw data for analysis, make sure is float for NaN values
        x_cl  = x_raw.astype(float)

        # flag 1005: data is not length N
        if (x_raw.size != N):
            flags.append(1005)

        # flag 1003: >1% are outside acceptable data range
        x_cl[np.where(x_raw < dataRng[0])] = np.nan
        x_cl[np.where(x_raw > dataRng[1])] = np.nan
        percOut = np.sum(np.isnan(x_cl))/float(N)
        if (percOut > 0.01):
              flags.append(1003)

        # flag 1007: sonic wind velocities are quantized
        if (is_quantized(x_raw) and any(x in datfield for \
                x in ['Sonic_u','Sonic_v','Sonic_w'])):
            flags.append(1007)

        # sort flags
        flags = sorted(flags)

        # remove spikes/detrend sonics, detrend everything else
        if (len(flags) == 0):
            if ('Sonic' in datfield):
                # remove spikes, linear detrend
                x_cl, n_spikes = cleantimeseries(t,x_cl)
                
                # flag 1008: too many spikes in sonic
                if (n_spikes >= n_spikes_max):
                    flags.append(1008)
                    
            else:
                x_cl = nandetrend(t,x_cl) + np.mean(x_cl)

        # save results in dictionary
        outdict          = {}
        outdict['time']  = t
        outdict['raw']   = x_raw
        outdict['clean'] = x_cl
        outdict['flags'] = flags

    elif (dataset == 'fluela'):

        # set data ranges, desired length of time series
        dataRng = dataRanges(dataset,datfield)
        dt, N   = 0.10, 6000
        t       = np.arange(N)*dt
                    
        # determine what kind of data was provided
        if (isinstance(data_in,dict)):                   # data structure given
            struc = data_in
        elif (isinstance(data_in,(int,float,tuple))):    # timestamp given
            fpath = time2fpath(dataset,data_in)
            struc = scio.loadmat(fpath)
        else:
            errStr = 'Incorrect data type ' + \
                '\"{}\" for input \"data_in\".'.format(type(data_in))
            raise TypeError(errStr)

        # initialize list of flags
        flags = []

        # try to load time series
        try:
            x_raw = np.squeeze(struc[datfield][0,0][0])            
            # flag 5003: data are all NaNs
            if (np.all(np.isnan(x_raw))):
                flags.append(5003)                
        except KeyError:
            # flag 1006: data are not in 20-Hz structure
            x_raw = np.empty(N)
            x_raw[:] = np.nan
            flags.append(1006)
                
        # copy raw data for analysis, make sure is float for NaN values
        x_cl  = x_raw.astype(float)

        # flag 1005: data is not length N
        if (x_raw.size != N):
            flags.append(1005)

        # flag 1003: >1% are outside acceptable data range
        x_cl[np.where(x_raw < dataRng[0])] = np.nan
        x_cl[np.where(x_raw > dataRng[1])] = np.nan
        percOut = np.sum(np.isnan(x_cl))/float(N)
        if (percOut > 0.01):
              flags.append(1003)

        # flag 1007: sonic wind velocities are quantized
        if (is_quantized(x_raw) and any(x in datfield for \
                x in ['Sonic_u','Sonic_v','Sonic_w'])):
            flags.append(1007)

        # sort flags
        flags = sorted(flags)

        # clean/detrend healthy data if no flags and it's a sonic time series
        if ((len(flags) == 0) and ('Sonic' in datfield)):
            x_cl = cleantimeseries(t,x_cl)

        # save results in dictionary
        outdict          = {}
        outdict['time']  = t
        outdict['raw']   = x_raw
        outdict['clean'] = x_cl
        outdict['flags'] = flags

    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return outdict


class tsin:
    """ A new class to track TurbSim input parameters of interest.
    """
    
    def __init__(self):
        self.ScaleIEC = [];
        self.n_z      = [];
        self.n_y      = [];
        self.dt       = [];
        self.t_anal   = [];
        self.t_use    = [];
        self.zhub     = [];
        self.turbc    = [];
        self.zref     = [];
        self.Vref     = [];


def readInput_v1(filename):
    """ Read specific parameters of interest from TurbSim .inp file and save in
        custom Python class tsin
        
        Args:
            filename (string): path to TS file to read
    """

    if ('.inp' in filename):

        # define line numbers for each parameter
        line_ScaleIEC = 15;
        line_n_z      = 18;
        line_n_y      = 19;
        line_dt       = 20;
        line_t_anal   = 21;
        line_t_use    = 22;
        line_zhub     = 23;
        line_turbc    = 32;
        line_zref     = 36;
        line_Vref     = 37;

        # initialize the class variable
        tsinput = tsin();

        # open file and extract parameters
        with open( filename, 'r' ) as f:
            idx = 1;
            for line in f:
                if ( idx == line_ScaleIEC ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.ScaleIEC = int(val);
                elif ( idx == line_n_z ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.n_z = int(val);
                elif ( idx == line_n_y ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.n_y = int(val);
                elif ( idx == line_dt ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.dt = float(val);
                elif ( idx == line_t_anal ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.t_anal = float(val);
                elif ( idx == line_t_use ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.t_use = float(val);
                elif ( idx == line_zhub ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.zhub = float(val);
                elif ( idx == line_turbc ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.turbc = val;
                elif ( idx == line_zref ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.zref = float(val);
                elif ( idx == line_Vref ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.Vref = float(val);
                idx += 1;

        return tsinput;

    elif ('.sum' in filename):

        print 'Cannot interpret .sum yet';
        return [];

    else:

        print 'That is neither a .sum or .inp';
        return [];


def readInput_v2(filename):
    """ Read specific parameters of interest from TurbSim .inp file and save in
        custom Python class tsin
        
        Args:
            filename (string): path to TS file to read
    """

    if ('.inp' in filename):

        # define line numbers for each parameter
        #  (correspond directly with numbers in text file)
        line_ScaleIEC = 16;
        line_n_z      = 19;
        line_n_y      = 20;
        line_dt       = 21;
        line_t_anal   = 22;
        line_t_use    = 23;
        line_zhub     = 24;
        line_turbc    = 34;
        line_zref     = 39;
        line_Vref     = 40;

        # initialize the class variable
        tsinput = tsin();

        # open file and extract parameters
        with open( filename, 'r' ) as f:
            idx = 1;
            for line in f:
                if ( idx == line_ScaleIEC ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.ScaleIEC = int(val);
                elif ( idx == line_n_z ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.n_z = int(val);
                elif ( idx == line_n_y ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.n_y = int(val);
                elif ( idx == line_dt ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.dt = float(val);
                elif ( idx == line_t_anal ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.t_anal = float(val);
##                elif ( idx == line_t_use ):
##                    val = line.lstrip(' ').split(' ')[0]\
##                          .replace('\"','');
##                    tsinput.t_use = float(val);
                elif ( idx == line_zhub ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.zhub = float(val);
                elif ( idx == line_turbc ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.turbc = val;
                elif ( idx == line_zref ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.zref = float(val);
                elif ( idx == line_Vref ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.Vref = float(val);
                idx += 1;

        return tsinput;

    elif ('.sum' in filename):

        print 'Cannot interpret .sum yet';
        return [];

    else:

        print 'That is neither a .sum or .inp';
        return [];

# %%============================================================================
# METADATA PROCESSING
# ==============================================================================


def metadataFields(dataset):
    """ Define list of fields to be stored in metadata table

        Args:
            dataset (string): flag for dataset

        Returns:
            fields (list): list of strings defining metadata
                columns
    """

    if (dataset == 'NREL'):
        fields = ['Record_Time','Processed_Time','Height','Wind_Speed_Cup', \
              'Wind_Direction','Precipitation', \
              'Mean_Wind_Speed', \
              'Sigma_u','Concentration_u','Location_u', \
              'Sigma_v','Concentration_v','Location_v', \
              'Sigma_w','Concentration_w','Location_w', \
              'up_wp','vp_wp','wp_Tp','up_vp','tau_u','tau_v','tau_w', \
              'MO_Length_interp','MO_Length_near']
    elif (dataset == 'fluela'):
        fields = ['Record_Time','Processed_Time','Height','Sonic_Cup', \
              'Sonic_Direction', 'Mean_Wind_Speed', \
              'Sigma_u','Concentration_u','Location_u', \
              'Sigma_v','Concentration_v','Location_v', \
              'Sigma_w','Concentration_w','Location_w', \
              'up_wp','vp_wp','wp_Tp','up_vp','tau_u','tau_v','tau_w', \
              'MO_Length']
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)
        
    return fields


def check_datfields(dataset,datfield):
    """ Define list of fields of all possible correct fieldnames for 
        error checking

        Args:
            dataset (string): flag for dataset
            datfield (string): dataset-specific fieldname

        Returns:
            fields (list): list of strings defining metadata
                columns
    """

    if (dataset == 'NREL'):
        fields = ['Sonic_u_15m','Sonic_u_30m','Sonic_u_50m','Sonic_u_76m',\
            'Sonic_u_100m','Sonic_u_131m','Sonic_v_15m','Sonic_v_30m',\
            'Sonic_v_50m','Sonic_v_76m','Sonic_v_100m','Sonic_v_131m',\
            'Sonic_w_15m','Sonic_w_30m','Sonic_w_50m','Sonic_w_76m',\
            'Sonic_w_100m','Sonic_w_131m','Sonic_Temp_rotated_15m',\
            'Sonic_Temp_rotated_30m','Sonic_Temp_rotated_50m',\
            'Sonic_Temp_rotated_76m','Sonic_Temp_rotated_100m',\
            'Sonic_Temp_rotated_131m','Air_Temp_3m','Air_Temp_26m',\
            'Air_Temp_88m','Cup_WS_10m','Cup_WS_26m','Cup_WS_80m',\
            'Cup_WS_88m','Cup_WS_134m','Vane_WD_10m','Vane_WD_26m',\
            'Vane_WD_88m','Vane_WD_134m','PRECIP_INTEN','Dewpt_Temp_3m',\
            'Dewpt_Temp_26m','Dewpt_Temp_88m','Dewpt_Temp_134m',\
            'Baro_Presr_3m']
    elif (dataset == 'fluela'):
        fields = ['Sonic_u_36m','Sonic_u_54m','Sonic_u_75m',\
                  'Sonic_v_36m','Sonic_v_54m','Sonic_v_75m',\
                  'Sonic_w_36m','Sonic_w_54m','Sonic_w_75m',\
                  'Sonic_Temp_rotated_36m','Sonic_Temp_rotated_54m',\
                  'Sonic_Temp_rotated_75m','Sonic_CupEqHorizSpeed_36m',\
                  'Sonic_CupEqHorizSpeed_54m','Sonic_CupEqHorizSpeed_75m',\
                  'Sonic_direction_36m','Sonic_direction_54m',\
                  'Sonic_direction_75m']
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)
        
    return datfield in fields


def datasetSpecs(dataset):
    """ Measurement specifics for dataset

        Args:
            dataset (string): flag for dataset

        Returns:
            n_t (int): number of time steps per record
            dt (float): time step
            heights (numpy array): measurement heights
    """
    import numpy as np

    if (dataset == 'NREL'):
        n_t     = 12000
        dt      = 0.05
        heights = np.array([15,30,50,76,100,131])
    elif (dataset == 'fluela'):
        n_t     = 6000
        dt      = 0.10
        heights = np.array([36,54,75])
    elif (dataset == 'PM06'):
        n_t     = 12000
        dt      = 0.05
        heights = []
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return (n_t,dt,heights)


def dataRanges(dataset,datfield):
    """ Acceptable range of time series

        Args:
            dataset (string): flag for dataset
            datfield (string): dataset-specific fieldname

        Returns:
            dataRng (list): [min,max] range of data
    """

    if (dataset == 'NREL'):

        # define data ranges
        if ('Sonic_u' in datfield):
            dataRng = [-35.05,35.05]
        elif ('Sonic_v' in datfield):
            dataRng = [-35.05,35.05]
        elif ('Sonic_w' in datfield):
            dataRng = [-29.95,29.95]
        elif ('Sonic_T' in datfield):
            dataRng = [-19.95,49.95]
        elif ('Air_Temp' in datfield):
            dataRng = [-50.,50.]
        elif ('Cup_WS' in datfield):
            dataRng = [0.,80.]
        elif ('Vane_WD' in datfield):
            dataRng = [0.,360.]
        elif ('PRECIP' in datfield):
            dataRng = [0.,4.]
        elif ('Dewpt' in datfield):
            dataRng = [-50.,50.]
        elif ('Presr' in datfield):
            dataRng = [740.,1000.]
        else:
              raise KeyError('Field {} not recognized.'.format(datfield))
    elif (dataset == 'fluela'):

        # define data ranges
        if ('Sonic_u' in datfield):
            dataRng = [-35.05,35.05]
        elif ('Sonic_v' in datfield):
            dataRng = [-35.05,35.05]
        elif ('Sonic_w' in datfield):
            dataRng = [-29.95,29.95]
        elif ('Sonic_T' in datfield):
            dataRng = [-19.95,49.95]
        elif ('Sonic_Cup' in datfield):
            dataRng = [  0.00,49.57]            # sqrt(2*(35.05**2))
        elif ('Sonic_direction' in datfield):
            dataRng = [ 0.00,360.00]            # angle from 0 to 360
        else:
              raise KeyError('Field {} not recognized.'.format(datfield))

    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return dataRng


def getBasedir(dataset):
    """ Get path to base directory and check if it exists

        Args:
            dataset (string): flag for dataset

        Returns:
            basedir (str): path to top level of data directory
    """
    import os
    import platform

    if (dataset == 'NREL'):
        if (platform.system() == 'Linux'):
            basedir = '/media/jrinker/JRinker SeaGate External/data/nrel-20Hz/'
        elif (platform.system() == 'Windows'):
#            basedir = 'G:\\data\\nrel-20Hz'
##            basedir = 'E:\\data\\nrel-20Hz'
            basedir = 'H:\\data\\nrel-20Hz'
        if not os.path.exists(basedir):
            errStr = 'Incorrect or unavailable base ' + \
                     'directory for dataset \"{}\".'.format(dataset)
            raise IOError(errStr)
    elif (dataset == 'fluela'):
        if (platform.system() == 'Linux'):
            basedir = '/media/jrinker/JRinker SeaGate External/data/fluela-high_freq/'
        elif (platform.system() == 'Windows'):
#            basedir = 'G:\\data\\fluela-high_freq'
##            basedir = 'E:\\data\\fluela-high_freq'
            basedir = 'H:\\data\\fluela-high_freq'
        if not os.path.exists(basedir):
            errStr = 'Incorrect or unavailable base ' + \
                     'directory for dataset \"{}\".'.format(dataset)
            raise IOError(errStr)

    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return basedir


def makemetadata(dataset):
    """ Construct metadata table
    """

    # get base directory, check it exists
    basedir = getBasedir(dataset)

    # process NREL dataset
    if (dataset == 'NREL'):
        metadata = makeNRELmetadata(basedir)

    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return metadata


def listmetadata(dataset,i,list_mats):
    """ Return list of parameters all heights for element i
        in list_mats
    """
    import numpy as np
    import sys, os
    import scipy.io as scio

    if (dataset in ['NREL','fluela']):
        heights  = datasetSpecs(dataset)[2]             # sampling heights
        basedir  = getBasedir(dataset)                  # base directory
        fpath    = os.path.join(basedir,list_mats[i])   # path to mat file
        n_fields = len(metadataFields(dataset))         # list wind parameters
        
        # try to load the structure, return arrays of NaNs if failed
        try:
            struc = scio.loadmat(fpath)
        except Exception as e:
            print('Cannot load {}'.format(fpath))
            print('  ' + str(e))
            parms    = np.empty(n_fields)
            parms[:] = np.nan
            return [parms for _ in range(heights.size)]

        # loop through sonic heights
        h_parms = []
        for height in heights:
            
            # calculate parameters
            row = struc2metadata(dataset,struc,height)

            # append to list of structure parametsr
            h_parms.append(row)
            
    else:
        errStr = 'Dataset {} not coded'.format(dataset)
        raise AttributeError(errStr)
                
    return h_parms


def struc2metadata(dataset,struc,height):
    """ Calculate metadata parameters from high-frequency .mat

        Args:
            struc (dictionary): 20-Hz structure
            height (int): measurement height

        Returns:
            parameters (numpy array): 1D array of metadata parameters
    """
    import numpy as np

    if (dataset in ['NREL','fluela']):

        # get list of metadata fieldnames
        md_fields = metadataFields(dataset)

        # intialize array of metadata parameters
        parameters = np.empty(len(md_fields))

        # calculate all fields
        outdict = calculatefield(dataset,struc,height)

        # assign fields to location in output array if dictionary is non-empty
        if (len(outdict) > 0):
            for i in range(len(md_fields)):
                field = md_fields[i]
                parameters[i] = outdict[field]
        else:
            parameters[:] = np.nan
            
    else:
        errStr = 'Dataset {} is not coded'.format(dataset)
        raise AttributeError(errStr)

    return parameters


def calculateKaimal(x,dt):
    """ Return the optimal Kaimal time scale tau found using a grid search. 
        Can be converted to length scale using L = tau * U.

        Args:
            x (numpy array): time history
            dt (float): time step
        
        Returns:
            tau (float): optimal Kaimal length scale
    """
    import numpy as np
    
    # if (2+)D array is fed in, halt with error
    if ((len(x.shape)>1) and (x.shape[0] != 1 and x.shape[1] != 1)):
        errStr = 'calculateKaimal only works on 1D arrays'
        raise ValueError(errStr)

    # squeeze to 1D, interpolate NaN values
    x_1D = np.squeeze(x)
    t    = np.arange(x.size)*dt
    notnan_idx = np.logical_not(np.isnan(x_1D))
    t_int = t[notnan_idx]
    x_int = x_1D[notnan_idx]
    x_1D  = np.interp(t,t_int,x_int)

    # grid search parameters
    n_m = 200;                                  # number of points in grid
    tau_l = -1;                                 # left tau coefficient
    tau_r = 4;                                  # right tau coefficient
    taugrid = np.logspace(tau_l,tau_r,n_m). \
        reshape(1,n_m)                          # tau search grid (row)

    # intermediate parameters
    n_t  = x_1D.size                            # no. of time steps
    T    = n_t*dt                               # total time
    n_f  = uniqueComponents(n_t)                # no. unique components
    df   = 1/T                                  # frequency resolution
    f    = df*np.arange(1,n_f). \
        reshape(n_f-1,1)                        # frequency vector (col)
    sig  = np.nanstd(x_1D)                      # std. deviation

    # get discrete spectrum from data
    X_dat  = np.fft.fft(x_1D). \
        reshape(n_t,1)/n_t                      # Fourier vector (row)
    Sk_dat = X2Sk(X_dat)[1:]* \
        np.ones((1,n_m))                        # data discrete PSD from f1 up

    # get n_f x n_m array of continuous spectra
    S_kaim = KaimalSpectrum(f,taugrid,sig)      # continuous Kaimal spectrum
    Sk_kaim = S_kaim * df                       # discrete Kaimal spectrum
    Sk_theo = np.empty(Sk_kaim.shape)
    for i in range(n_m):
        alpha = spectralScale(Sk_kaim[:,i],sig,n_t)
##        alpha = 1
        Sk_theo[:,i] = (alpha**2)*Sk_kaim[:,i]
    
    J = np.sum(np.power(Sk_dat-Sk_theo, \
        2),axis=0)                              # array of squared errors

    i_min = np.argmin(J)                        # index of min error

    tau = taugrid[0,i_min]                      # tau_minError
    
    return tau


def interpolationHeights(dataset,ht,field):
    """ Get upper and lower measurement heights for interpolation
    
        Args:
            dataset (string): flag to indicate which dataset to analyze
            struc (numpy structure): loaded from matlab
            ht (float): measurement height to interpolate
            field (string): field of interest for interpolation

        Returns:
            interp_heights (tuple): (lowerHt,upperHt)
            
    """
    import numpy as np
    
    if (dataset == 'NREL'):

        if (field == 'Wind_Speed_Cup'):
            msmnt_hts = np.array([10,26,80,88,134])

        elif (field == 'Temperature'):
            msmnt_hts = np.array([3,26,88])
            
        elif (field == 'Wind_Direction'):
            msmnt_hts = np.array([10,26,88,134])
            
        elif (field == 'Dewpt_Temp'):
            msmnt_hts = np.array([3,26,88,134])

        else:
            errStr = 'Field \"{}\" is not coded '.format(field) + \
                'for dataset \"NREL\".'
            raise AttributeError(errStr)    
                        
        # force height to be in range
        ht = min(msmnt_hts[-1],ht)
        ht = max(msmnt_hts[0],ht)
        
        # find upper and lower heights
        upperHt = msmnt_hts[np.where(ht <= msmnt_hts)[0][0]]
        lowerHt = msmnt_hts[np.where(ht <= msmnt_hts)[0][0] - 1]
        
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)
        
    return (lowerHt,upperHt)
    
    
def interpolateparameter(dataset,ht,lo_val,hi_val,field):
    """ Interpolate parameter value
    """
    import numpy as np
    
    if (dataset == 'NREL'):
        
        # get heights for that measurement
        lo_ht, hi_ht = interpolationHeights(dataset,ht,field)
        
        # linearly interpolate linear variables
        if (field == 'Wind_Speed_Cup'):
            msmnt_hts = np.array([10,26,80,88,134])
            ht = min(msmnt_hts[-1],ht)          # force height to be in 
            ht = max(msmnt_hts[0],ht)           #   measurement range
            val = (hi_val - lo_val) / float((hi_ht - lo_ht)) * \
                (ht - lo_ht) + lo_val
                
        # linearly interpolate linear variables
        elif (field == 'Temperature'):
            msmnt_hts = np.array([3,26,88])
            ht = min(msmnt_hts[-1],ht)          # force height to be in 
            ht = max(msmnt_hts[0],ht)           #   measurement range
            val = (hi_val - lo_val) / float((hi_ht - lo_ht)) * \
                (ht - lo_ht) + lo_val
            
        # linearly interpolate linear variables
        elif (field == 'Dewpt_Temp'):
            msmnt_hts = np.array([3,26,88,134])
            ht = min(msmnt_hts[-1],ht)          # force height to be in 
            ht = max(msmnt_hts[0],ht)           #   measurement range
            val = (hi_val - lo_val) / float((hi_ht - lo_ht)) * \
                (ht - lo_ht) + lo_val
            
        # special interpolate wrapping variables
        elif (field == 'Wind_Direction'):
            msmnt_hts = np.array([10,26,88,134])
            ht = min(msmnt_hts[-1],ht)          # force height to be in 
            ht = max(msmnt_hts[0],ht)           #   measurement range
            dtheta = np.angle(np.exp(1j*hi_val*np.pi/180.)/\
                np.exp(1j*lo_val*np.pi/180.),deg=1)
            val = (dtheta / (hi_ht - lo_ht) * (ht - lo_ht) + lo_val) % 360.
                        
        else:
            raise KeyError('Unknown field \"{}\" for '.format(field) + \
                'dataset \"NREL\"')
                        
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)
        
    return val


def calculatefield(dataset,struc20,ht):
    """ Save atmophseric parameters in output dictionary

    """
    import calendar, time
    import numpy as np
    
    # meteorological constants
    g, R, kappa = 9.81, 287, 0.41

    if (dataset == 'NREL'):

        dt, N = 0.05, 12000
        refht = 3
        
        # get interpolation/closest heights
        T_hts = interpolationHeights(dataset,ht,'Temperature')
        loht_T, hiht_T = T_hts
        clht_T = T_hts[np.abs(np.array(T_hts)-ht).argmin()]
        loht_WS, hiht_WS = interpolationHeights(dataset,ht,'Wind_Speed_Cup')
        loht_WD, hiht_WD = interpolationHeights(dataset,ht,'Wind_Direction')
        loht_DP, hiht_DP = interpolationHeights(dataset,ht,'Dewpt_Temp')
        
        # load all necessary time series, incrementally checking flags
        clean  = 1
        fields = ['Sonic_u','Sonic_v','Sonic_w','Sonic_T','Temperature',\
            'Temperature','Temperature','Wind_Speed_Cup','Wind_Speed_Cup',\
            'Wind_Direction','Wind_Direction','Precipitation',\
            'Dewpt_Temp','Dewpt_Temp','Dewpt_Temp','Pressure',\
            'Temperature']
        ts_hts = [ht, ht, ht, ht, loht_T, hiht_T, clht_T, loht_WS, hiht_WS,\
                loht_WD, hiht_WD, ht, loht_DP, hiht_DP, refht, refht,\
                refht]
        time_series = np.empty((len(fields),N))
        for i in range(len(fields)):
            
            # load time series for that field, measurement height
            field = fields[i]
            ts_ht = ts_hts[i]
            outdict  = loadtimeseries(dataset,field,ts_ht,struc20)

            # save time series if it is clean
            if (len(outdict['flags']) == 0):
                time_series[i,:] = outdict['clean']
            else:
                clean = 0
                
        # if all time series are clean
        if clean:
            
            # put all time series in variables
            t     = np.arange(N)*dt                 # time vector
            u     = time_series[0,:]                # Sonic_u
            v     = time_series[1,:]                # Sonic_v
            w     = time_series[2,:]                # Sonic_w
            T_s   = time_series[3,:]                # Sonic_T
            T_lo  = time_series[4,:]                # lower ht temp
            T_hi  = time_series[5,:]                # upper ht temp
            T_cl  = time_series[6,:]                # closest ht temp
            WS_lo = time_series[7,:]                # lower ht cup WS
            WS_hi = time_series[8,:]                # upper ht cup WS
            WD_lo = time_series[9,:]                # lower ht wind dir
            WD_hi = time_series[10,:]               # upper ht wind dir
            P     = time_series[11,:]               # precipitation
            DP_lo = time_series[12,:]               # lower ht dewpoint temp
            DP_hi = time_series[13,:]               # upper ht dewpoint temp
            DP_0  = time_series[14,:]               # ref ht dewpoint temp 
            P_0   = time_series[15,:]               # ref ht pressure
            T_0   = time_series[16,:]               # ref ht temperature
            
            # get mean values for later calculations
            fname     = struc20['tower'][0,0][23][0,0][0][0,0][8][0]
            rec_time  = NRELfname2time(fname)
            up        = nandetrend(t,u)
            vp        = nandetrend(t,v)
            wp        = nandetrend(t,w)
            Tp        = nandetrend(t,T_s)
            upwp_bar  = np.nanmean(up*wp)
            upvp_bar  = np.nanmean(up*vp)
            vpwp_bar  = np.nanmean(vp*wp)
            wpTp_bar  = np.nanmean(wp*Tp)
            WSbar_lo  = np.nanmean(WS_lo)
            WSbar_hi  = np.nanmean(WS_hi)
            WDbar_lo  = np.angle(np.nanmean(np.exp(1j*WD_lo*np.pi/180.)),deg=1)
            WDbar_hi  = np.angle(np.nanmean(np.exp(1j*WD_hi*np.pi/180.)),deg=1)
            DPbar_lo  = np.nanmean(DP_lo)
            DPbar_hi  = np.nanmean(DP_hi)
            Pbar_0    = np.nanmean(P_0)
            Tbar_0    = np.nanmean(T_0)
            DPbar_0   = np.nanmean(DP_0)

            # calculate derived intermediate values
            ustar     = ( (upwp_bar)**2 + (vpwp_bar)**2 ) ** 0.25
            rhou,muu  = signalPhaseCoherence(up)
            rhov,muv  = signalPhaseCoherence(vp)
            rhow,muw  = signalPhaseCoherence(wp)
            Tbar_lo_K = C2K(np.nanmean(T_lo))
            Tbar_hi_K = C2K(np.nanmean(T_hi))
            Tbar_cl_K = C2K(np.nanmean(T_cl))
            Tbar_in_K = interpolateparameter(dataset,ht,Tbar_lo_K,Tbar_hi_K,\
                                            'Temperature')
            WSbar_in  = interpolateparameter(dataset,ht,WSbar_lo,WSbar_hi,\
                                            'Wind_Speed_Cup')
            WDbar_in  = interpolateparameter(dataset,ht,WDbar_lo,WDbar_hi,\
                                            'Wind_Direction')
            Tbar_0_K  = C2K(Tbar_0)
            if (DPbar_0 > 0): A, B = 7.5, 237.3
            else:            A, B = 9.5, 265.5
            e0        = 6.11 * 10 ** ((DPbar_0*A)/(DPbar_0 + B))
            q0        = e0 / Pbar_0
            Tv0       = (Tbar_0_K)*(1 + 0.61*q0)
            dPdz      = - (g * Pbar_0) / (R * Tv0)
            DPbar_z   = interpolateparameter(dataset,ht,DPbar_lo,
                                             DPbar_hi,'Dewpt_Temp')
            Tbar_z_K  = interpolateparameter(dataset,ht,Tbar_lo_K,
                                             Tbar_hi_K,'Temperature')

            if (DPbar_z > 0): A, B = 7.5, 237.3
            else:            A, B = 9.5, 265.5
            ez        = 6.11 * 10 ** ((DPbar_z*A)/(DPbar_z + B))
            Pz        = Pbar_0 + (ht - refht)*dPdz
            qz        = ez / Pz
            Tvbar_z_K = (Tbar_z_K)*(1 + 0.61*qz)

            # initialize output dictionary
            outdict = {}

            # save values
            outdict['Record_Time']     = rec_time
            outdict['Processed_Time']  = calendar.timegm(time.gmtime())   
            outdict['Height']          = ht
            outdict['Wind_Speed_Cup']  = WSbar_in
            outdict['Wind_Direction']  = WDbar_in
            outdict['Precipitation']   = np.nanmean(P)
            outdict['Mean_Wind_Speed'] = np.nanmean(u)
            outdict['Sigma_u']         = np.nanstd(up)
            outdict['Concentration_u'] = rhou
            outdict['Location_u']      = muu
            outdict['Sigma_v']         = np.nanstd(vp)
            outdict['Concentration_v'] = rhov
            outdict['Location_v']      = muv
            outdict['Sigma_w']         = np.nanstd(wp)
            outdict['Concentration_w'] = rhow
            outdict['Location_w']      = muw
            outdict['up_wp']           = upwp_bar          
            outdict['vp_wp']           = vpwp_bar
            outdict['wp_Tp']           = wpTp_bar
            outdict['up_vp']           = upvp_bar
            outdict['tau_u']           = calculateKaimal(up + \
                                                np.nanmean(u),dt)
            outdict['tau_v']           = calculateKaimal(vp,dt)
            outdict['tau_w']           = calculateKaimal(wp,dt)
            outdict['MO_Length_interp'] = -(Tbar_in_K * ustar**3) \
                                         /(kappa * g * wpTp_bar)
            outdict['MO_Length_near']  = - (Tbar_cl_K * ustar**3)/ \
                                         (kappa * g * wpTp_bar)
            outdict['MO_Length_virt']  = - (Tvbar_z_K * ustar**3)/ \
                                         (kappa * g * wpTp_bar)
            
        else:
            outdict = {}

    elif (dataset == 'fluela'):

        dt, N = 0.10, 6000
        
        # no need to interpolate -- only have sonic measurements
        
        # load all necessary time series, incrementally checking flags
        clean  = 1
        fields = ['Sonic_u','Sonic_v','Sonic_w','Sonic_T',\
            'Sonic_Cup','Sonic_Direction']
        ts_hts = [ht, ht, ht, ht, ht, ht]
        time_series = np.empty((len(fields),N))
        for i in range(len(fields)):
            
            # load time series for that field, measurement height
            field = fields[i]
            ts_ht = ts_hts[i]
            outdict  = loadtimeseries(dataset,field,ts_ht,struc20)
            
            # save time series if it is clean
            if (len(outdict['flags']) == 0):
                time_series[i,:] = outdict['clean']
            else:
                clean = 0
                
        # if all time series are clean
        if clean:
            
            # put all time series in variables
            t     = np.arange(N)*dt                 # time vector
            u     = time_series[0,:]                # Sonic_u
            v     = time_series[1,:]                # Sonic_v
            w     = time_series[2,:]                # Sonic_w
            T_s   = time_series[3,:]                # Sonic_T
            WS    = time_series[4,:]                # Sonic Cup WS
            WD    = time_series[5,:]                # Sonic Wind Direction

            
            # get variables necessary for later calculations
            fname     = struc20['tower'][0,0][23][0,0][0][0,0][3][0]
            rec_time  = NRELfname2time(fname)
            T_s_K     = C2K(T_s)
            up        = nandetrend(t,u)
            vp        = nandetrend(t,v)
            wp        = nandetrend(t,w)
            Tp        = nandetrend(t,T_s_K)
            upwp_bar  = np.nanmean(up*wp)
            upvp_bar  = np.nanmean(up*vp)
            vpwp_bar  = np.nanmean(vp*wp)
            wpTp_bar  = np.nanmean(wp*Tp)
            ustar     = ( (upwp_bar)**2 + (vpwp_bar)**2 ) ** 0.25
            rhou,muu  = signalPhaseCoherence(up)
            rhov,muv  = signalPhaseCoherence(vp)
            rhow,muw  = signalPhaseCoherence(wp)
            Tbar_K    = np.nanmean(T_s_K)
            WSbar     = np.nanmean(WS)
            WDbar     = np.angle(np.nanmean(np.exp(1j*WD*np.pi/180.)),deg=1)
    
            # initialize output dictionary
            outdict = {}
    
            # save values
            outdict['Record_Time']     = rec_time
            outdict['Processed_Time']  = calendar.timegm(time.gmtime())   
            outdict['Height']          = ht
            outdict['Sonic_Cup']       = WSbar
            outdict['Sonic_Direction'] = WDbar
            outdict['Mean_Wind_Speed'] = np.nanmean(u)
            outdict['Sigma_u']         = np.nanstd(up)
            outdict['Concentration_u'] = rhou
            outdict['Location_u']      = muu
            outdict['Sigma_v']         = np.nanstd(vp)
            outdict['Concentration_v'] = rhov
            outdict['Location_v']      = muv
            outdict['Sigma_w']         = np.nanstd(wp)
            outdict['Concentration_w'] = rhow
            outdict['Location_w']      = muw
            outdict['up_wp']           = upwp_bar          
            outdict['vp_wp']           = vpwp_bar
            outdict['wp_Tp']           = wpTp_bar
            outdict['up_vp']           = upvp_bar
            outdict['tau_u']           = calculateKaimal(up + \
                                                np.nanmean(u),dt)
            outdict['tau_v']           = calculateKaimal(vp,dt)
            outdict['tau_w']           = calculateKaimal(wp,dt)
            outdict['MO_Length']       = -(Tbar_K * ustar**3) \
                                         /(0.41 * 9.81 * wpTp_bar)
            
        else:
            outdict = {}
        
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return outdict


def list_matfiles(basedir,save=0):
    """ List all .mat files in subdirectories of base directory.
        Save in file in top level of base directory if desired.

        Args:
            basedir (string): path to top level of base directory
            save (boolean): opt, flag to save list

        Returns:
            list_mats (list): list of paths to mat files
    """
    import sys, os, datetime, json

    list_mats = []
    print('Processing files...')
    sys.stdout.flush()
    for root, dirs, files in os.walk(basedir,topdown=False):
        for fname in files:
            if fname.endswith('.mat'):
                fpath     = os.path.join(root,fname)
                loc_fpath = fpath[len(basedir):].lstrip(os.path.sep)
                list_mats.append(loc_fpath)
    print('Completed')

    if save:
        fname = 'listmats_' + datetime.date.today(). \
                isoformat() + '.txt'
        fpath = os.path.join(basedir,fname)
        with open(fpath,'w') as f:
            json.dump(list_mats,f)
        print('Data saved to ' + fpath)

    return list_mats

# %%============================================================================
# METADATA ANALYSIS
# ==============================================================================


def metadataFpath(dataset):
    """ Filepath to metadata table
    
        Args:
            dataset (string): toggle for dataset choice
            
        Returns:
            fpath (string): filepath to metadata file
    """
    import os
    
    base  =  'C:\\Users\\jrinker\\Dropbox\\research\\processed_data'
    fpath = os.path.join(base,dataset + '-metadata.mat')
        
    return fpath


def screenmetadata(fields,metadata,dataset):
    """ Screen the metadata for data quality
    
        Args:
            fields (list): list of fields associated with metadata cols
            metadata(numpy array): array of metadata
            dataset (string): toggle for dataset choice
            
        Returns:
            cleandata (numpy array): numpy array with screened data
    """
    import numpy as np
    
    if (dataset == 'NREL'):
        CSLim  = 3                          # lower cup speed limit
        dir1   = 240                        # CCW edge for direction range
        dir2   = 315                        # CW edge for direction range
        preLim = 2.7                        # lower precipitation limit
        
        # column indices for each value
        CScol  = fields.index('Wind_Speed_Cup')
        dirCol = fields.index('Wind_Direction')
        preCol = fields.index('Precipitation')
        
        # filter out the rows with NaN values
        metadata = metadata[np.logical_not( \
            np.any(np.isnan(metadata),axis=1)),:]
        
        # screen remaining data
        cleandata = metadata[np.where(metadata[:,CScol] > CSLim)]
        cleandata = cleandata[np.where(cleandata[:,dirCol] >= dir1)]
        cleandata = cleandata[np.where(cleandata[:,dirCol] <= dir2)]
        cleandata = cleandata[np.where(cleandata[:,preCol] >= preLim)]

    elif (dataset == 'fluela'):
        CSLim = 3                           # lower cup speed limit
        dir1  = -180                        # CCW edge for direction range
        dir2  = 0                           # CW edge for direction range
        T1    = timetup2flt((2010,1,1,0,0)) # start time
        T2    = timetup2flt((2010,3,23,0,0))# end time
        
        # column indices for each value
        CScol   = fields.index('Sonic_Cup')
        dirCol  = fields.index('Sonic_Direction')
        timeCol = fields.index('Record_Time')
        
        # filter out the rows with NaN values
        metadata = metadata[np.logical_not( \
            np.any(np.isnan(metadata),axis=1)),:]
        
        # screen remaining data
        cleandata = metadata[np.where(metadata[:,CScol] > CSLim)]
        cleandata = cleandata[np.where(cleandata[:,dirCol] >= dir1)]
        cleandata = cleandata[np.where(cleandata[:,dirCol] <= dir2)]
        cleandata = cleandata[np.where(cleandata[:,timeCol] >= T1)]
        cleandata = cleandata[np.where(cleandata[:,timeCol] <= T2)]
        
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)
        
    return cleandata


def compositeCDF(x,dist_name,p_main,x_T=float("inf"),p_GP=(0.1,0,1)):
    """ Cumulative distribution function of single/composite distribution

        Args:
            x (numpy array): values at which to evaluate CDF
            dist_name (string): distribution type
            p_main (tuple): main distribution parameters
            x_T (float): optional, threshold value
            p_GP (tuple): optional, GP distribution parameters

        Returns:
            F (numpy array): composite CDF values at x
    """
    import numpy as np
    import scipy.stats
    
    # initialize cdf
    F_comp = np.empty(x.shape)
    
    # initialize distributions
    dist_main = getattr(scipy.stats, dist_name)
    dist_GP   = getattr(scipy.stats, 'genpareto')
    
    # define CDF functions to improve code readability
    F_main = lambda x: dist_main.cdf(x, *p_main[:-2], \
            loc=p_main[-2], scale=p_main[-1])
    F_GP   = lambda x: dist_GP.cdf(x, *p_GP[:-2], \
            loc=p_GP[-2], scale=p_GP[-1])
    
    # get indices for main and GP distributions
    i_main, i_GP = np.where(x <= x_T), np.where(x > x_T)
    
    # set distribution values
    F_comp[i_main] = F_main(x[i_main])
    F_comp[i_GP]   = (1 - F_main(x_T)) * F_GP(x[i_GP]) + F_main(x_T)
    
    return F_comp


def inversecompositeCDF(Q,dist_name,p_main,x_T=float("inf"),p_GP=(0.1,0,1)):
    """ Inverse cumulative distribution function of single/composite
        distribution

        Args:
            Q (numpy array): values at which to evaluate CDF
            dist_name (string): distribution type
            p_main (tuple): main distribution parameters
            x_T (float): optional, threshold value
            p_GP (tuple): optional, GP distribution parameters

        Returns:
            x (numpy array): composite CDF values at x
    """
    import numpy as np
    import scipy.stats
    
    # initialize array for output
    x_comp = np.empty(Q.shape)
    
    # initialize distributions
    dist_main = getattr(scipy.stats, dist_name)
    dist_GP   = getattr(scipy.stats, 'genpareto')
    
    # define CDF functions to improve code readability
    G_main = lambda Q: dist_main.ppf(Q, *p_main[:-2], \
            loc=p_main[-2], scale=p_main[-1])
    G_GP   = lambda Q: dist_GP.ppf(Q, *p_GP[:-2], \
            loc=p_GP[-2], scale=p_GP[-1])

    # calculate threshold quantile
    Q_T = dist_main.cdf(x_T, *p_main[:-2], \
            loc=p_main[-2], scale=p_main[-1])
    
    # get indices for main and GP distributions
    i_main, i_GP = np.where(Q <= Q_T), np.where(Q > Q_T)
    
    # set distribution values
    x_comp[i_main] = G_main(Q[i_main])
    x_comp[i_GP]   = G_GP((Q[i_GP]-Q_T)/(1-Q_T))
    
    return x_comp

    
def compositePDF(x,dist_name,p_main,x_T=float("inf"),p_GP=(0.1,0,1)):
    """ Probability density function of single/composite distribution

        Args:
            x (numpy array): values at which to evaluate CDF
            dist_name (string): distribution type
            p_main (tuple): main distribution parameters
            x_T (float): optional, threshold value
            p_GP (tuple): optional, GP distribution parameters

        Returns:
            f (numpy array): composite PDF values at x
    """
    import numpy as np
    import scipy.stats
    
    # initialize pdf
    f_comp = np.empty(x.shape)
    
    # initialize distributions
    dist_main = getattr(scipy.stats, dist_name)
    dist_GP   = getattr(scipy.stats, 'genpareto')
    
    # define PDF functions to improve code readability
    f_main = lambda x: dist_main.pdf(x, *p_main[:-2], \
            loc=p_main[-2], scale=p_main[-1])
    f_GP   = lambda x: dist_GP.pdf(x, *p_GP[:-2], \
            loc=p_GP[-2], scale=p_GP[-1])
    
    # get indices for main and GP distributions
    i_main, i_GP = np.where(x <= x_T), np.where(x > x_T)
    
    # set distribution values
    f_comp[i_main] = f_main(x[i_main])
    f_comp[i_GP]   = f_GP(x[i_GP])
    
    return f_comp


def compositeNSAE(x,dist_name,p_main,x_T=float("inf"),p_GP=(0.1,0,1)):
    """ Normalized sum of absolute error for single/composite distribution

        Args:
            x (numpy array): values at which to evaluate CDF
            dist_name (string): distribution type
            p_main (tuple): main distribution parameters
            x_T (float): optional, threshold value
            p_GP (tuple): optional, GP distribution parameters

        Returns:
            NSAE (float): normalized sum absolute error
    """
    import numpy as np
    
    # get threshold value
    N = x.size
    x = np.sort(x)
    
    # get CDFs
    F_emp = np.arange(1,N+1)/(N+1.)
    F_fit = compositeCDF(x,dist_name,p_main,x_T,p_GP)
    
    # calculate NSAE
    NSAE  = np.mean(np.abs(F_fit-F_emp))
    
    return NSAE


def fitcompositeparameters(x,dist_name,x_T=float('inf')):
    """ Optimize single/composite parameters by minimizing NSAE, given
        threshold value (if x_T omitted, fits single distribution)

        Args:
            x (numpy array): data
            x_T (float): threshold between distributions
            dist_name (string): main distribution

        Returns:
            p_main_opt (tuple): main distribution parameters
            p_GP_opt (tuple): GP distribution parameters
    """
    import scipy.stats
    from scipy.optimize import minimize
    import numpy as np

    # initialize CDFs
    dist_main = getattr(scipy.stats, dist_name)
    dist_GP   = getattr(scipy.stats, 'genpareto')

    # get initial guesses using MLE
    p_main0 = dist_main.fit(x)
    p_GP0   = dist_GP.fit(x[np.where(x > x_T)], loc=x_T)
    p_comp0 = p_main0 + (p_GP0[0],) + (p_GP0[-1],)

    # get parameter bounds
    bnds_main = parameterbounds(dist_name)
    bnds_GP   = parameterbounds('genpareto')
    bnds      = bnds_main + (bnds_GP[0],) + (bnds_GP[-1],)

    # define error function
    fun = lambda p_comp: compositeNSAE(x,dist_name,p_comp[:-2],x_T, \
        (p_comp[-2],)+(x_T,)+(p_comp[-1],))

    # perform bounded optimization
    res = minimize(fun,p_comp0,bounds=bnds)

    # convert numpy array to tuple
    param_opt = tuple(np.ndarray.tolist(res.x))

    # separate main distribution and GP distribution parameters
    p_main_opt = param_opt[:-2]
    p_GP_opt   = (param_opt[-2],) + (x_T,) + (param_opt[-1],)

    return (p_main_opt,p_GP_opt)


def parameterbounds(dist_name):
    """ Tuple of bounds on parameters for given distribution name

        Args:
            dist_name (string): Python name of distribution

        Returns:
            bnds (tuple): (min,max) pairs for each parameter
    """

    if (dist_name == 'anglit'):
        bnds = ((None,None),(0,None))
    elif (dist_name == 'arcsine'):
        bnds = ((None,None),(0,None))
    elif (dist_name == 'beta'):
        bnds = ((0,None),(0,None),(None,None),(0,None))
    elif (dist_name == 'chi2'):
        bnds = ((0,None),(None,None),(0,None))
    elif (dist_name == 'cosine'):
        bnds = ((None,None),(0,None))
    elif (dist_name == 'expon'):
        bnds = ((None,None),(0,None))
    elif (dist_name == 'exponweib'):
        bnds = ((0,None),(0,None),(None,None),(0,None))
    elif (dist_name == 'genextreme'):
        bnds = ((None,None),(None,None),(0,None))
    elif (dist_name == 'gengamma'):
        bnds = ((0,None),(None,None),(None,None),(0,None))
    elif (dist_name == 'genpareto'):
        bnds = ((None,None),(None,None),(0,None))
    elif (dist_name == 'lognorm'):
        bnds = ((0,None),(None,None),(0,None))
    elif (dist_name == 'vonmises'):
        bnds = ((0,None),(None,None),(0,None))
    elif (dist_name == 'wrapcauchy'):
        bnds = ((0,None),(None,None),(0,None))
    else:
        print('***THAT DISTRIBUTION IS NOT CODED***')
        bnds = ()

    return bnds


def fitcompositedistribution(dataset,iP,x):
    """ Optimize NSAE, returning distribution type, distribution parameters,
        threshold value, and Generalized Pareto distribution

        Args:
            dataset (string): flag to indicate source of data
            iP (integer): parameter index (0=U,1=s_u,2=L,3=rho)
            x (numpy array): data

        Returns:
        
    """
    import numpy as np

    if (dataset == 'NREL'):

        # probability distribution candidates
        half_cands = ['lognorm','genextreme',\
        'chi2','gengamma','exponweib','expon']      # U, sigma_u, L
        fine_cands = ['anglit','beta', \
            'genextreme','gengamma','lognorm']      # rho

        # quantile threshold values to search
        Q_Ts = [0.80,0.85,0.90,0.95,1.00]

        # set distribution candidates
        if (iP < 3):
            dist_cands = half_cands
        elif (iP == 3):
            dist_cands = fine_cands

        # initialize error matrix
        NSAEs = np.empty((len(dist_cands),len(Q_Ts)))

        # initialize parameter list
        d_parms = []

        # loop through distribution candidates
        for iD in range(len(dist_cands)):
            
            # isolate distribution
            dist_name = dist_cands[iD]
            print('Processing {}'.format(dist_name))

            # initialize parameter list
            Q_parms = []

            # loop through threshold values
            for iQ in range(len(Q_Ts)):

                # isolate quantile
                Q_T  = Q_Ts[iQ]
                print('  ...quantile {}'.format(Q_T))

                # calculate threshold value
                x, N = np.sort(x), x.size
                if (Q_T == 1.0): x_T = float('inf')
                else:            x_T  = x[N*Q_T]

                # optimize distribution parameters
                p_main, p_GP = fitcompositeparameters( \
                    x,dist_name,x_T)

                # calculate corresponding NSAE
                NSAE = compositeNSAE(x,dist_name,p_main, \
                                     x_T,p_GP)

                # save parameters and NSAE
                fit_parms    = (dist_name,p_main,x_T,p_GP,NSAE)
                NSAEs[iD,iQ] = NSAE

                Q_parms.append(fit_parms)

            d_parms.append(Q_parms)

    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    # return values corresponding to minimum NSAE
    iD_min, iQ_min = np.where(NSAEs == NSAEs.min())
    fit_parms_min = d_parms[iD_min][iQ_min]

    print(NSAEs)
    print(iD_min,iQ_min)

    return fit_parms_min


# ==============================================================================
# SIMULATION
# ==============================================================================


def wrappedCauchySample(n_t,n_m,rho,mu):
    """ Random numbers from a Wrapped Cauchy distribution (n_t x n_m) with 
        location parameter mu and concentration parameter rho.
        
        Reference: Statistical Analysis of Circular Data, Fisher, Sec. 3.3.4
        Modified to correctly implement arccos.
        
        Args:
            n_t (int): sample size
            n_m (int): number of samples
            rho (float): concentration parameter
            mu (float): location parameter
    
        Returns:
            theta (numpy array): sample of angles
    """
    import numpy as np
    
    if ((rho < 0) or (rho>1)):
        print 'rho must be between 0 and 1'
        return []
    
    U = np.random.rand(n_t,n_m)                     # uniform random numbers
    
    V = np.cos(2.*np.pi*U)                          # See lines before Eq. 3.28
    c = 2.*rho/(1+(rho**2))                         # See lines before Eq. 3.28
    
    B = 2*np.round(np.random.rand(n_t,n_m)) - 1     # boolean RV
                
    theta = wrap(np.multiply(B, (np.arccos(np.divide(V+c, \
        1+c*V)))) + mu)                             # sample of angles
    
    return theta


def generateKaimal1D(n_t,n_m,dt,U,sig,tau,rho,mu):
    """ Single-point turbulent time history with the Kaimal spectrum and the 
        given parameters.

        Args:
            n_t (int): number time steps
            n_m (int): number realizations
            dt (float): time step size
            U (float): mean wind speed
            sig (float): wind standard deviation
            tau (float): Kaimal integral time scale
            rho (float): concentration parameter
            mu (float): location parameter
        
        Returns:
            t (numpy array): time
            x (numpy array): [n_t x n_m] array of turbulent wind records
    """
    import numpy as np

    scale = 1;    # flag to scale for time discretizaion

    # intermediate parameters
    T    = n_t*dt;                                  # total time
    t    = np.arange(0,n_t)*dt;                     # time vector
    n_f  = uniqueComponents(n_t);                   # no. unique components
    df   = 1./T;                                    # frequency resolution
    f    = df*np.arange(0,n_f);                     # frequency array

    if (sig<0):
        print 'Turbulence cannot be less than zero.';
        return []
        
    elif (sig==0):
        x = np.ones((n_t,n_m))*U;
        return [t,x]
        
    else:
        S = KaimalSpectrum(f,tau,sig);              # Kaimal spectrum
        Sk = S*df;                                  # discrete spectrum
        Xmag = Sk2Xuniq(Sk).reshape(n_f,1)          # S[k] -> |X[k]|
        
        if scale:                                   # spectral scaling
            alpha = spectralScale(Sk,sig,n_t);      #   for discretizaion
            Xmag = alpha*Xmag

        dtheta = wrappedCauchySample(n_f-1,\
            n_m,rho,mu)                             # phase differences
        phases = np.cumsum(dtheta,axis=0)           # phases from f1 to end
        phases = np.append(np.zeros((1,n_m)), \
            phases, axis=0)                         # phases from f0 to end
                
        Xpha = np.exp(1j*phases)                    # phases -> phasors
        
        Xuniq = np.multiply(Xmag,Xpha)              # X = mags*phasors
        Xuniq[0][:] = U;                            # real signal, correct mean
        X = Xuniq2X(Xuniq,n_t)                      # full Fourier vector
        
        x = np.real(np.fft.ifft(X,axis=0) \
            *n_t).reshape(n_t,n_m)                  # time series
                
        return (t,x)


def IEC_VelProfile(z,zhub,Vhub):
    """ IEC velocity profile for values in array z with hub height zhub and hub
        velocity Vhub.
        
        Args:
            z (numpy array): array of heights in meters for velocity calcs
            zhub (float): hub-height of turbine in meters
            Vhub (float): hub-height mean velocity in m/s
        
        Returns:
            V (numpy array): array of mean wind speeds at heights in z
    """
    import numpy as np

    alpha = 0.2;

    V = Vhub * np.power( \
        z/zhub, alpha );

    return V


def IEC_Iref(turbc):
    """ Turbulence intensity for turbulence class turbc
    
        Args:
            turbc (string): turbulence class from IEC 61400-1 Ed. 3
        
        Returns:
            Iref (float): reference turbulence intensity
    """

    if ('A' in turbc):                  # class A
        Iref = 0.16;
    elif ('B' in turbc):                # class B
        Iref = 0.14;
    elif ('C' in turbc):                # class C
        Iref = 0.12;
    else:                               # user's choice
        Iref = float(turbc);

    return Iref


def IEC_Sigma1(Iref,Vhub):
    """ Longitudinal standard deviation
    
        Args:
            Iref (float): reference hub-height turbulence intensity
            
        Returns:
            sigma1 (float): longitudinal standard deviation
    """
    
    sigma1 = Iref*(0.75 * Vhub + 5.6);
    
    return sigma1


def IEC_TiProfile(z,zhub,Vhub,turbc):
    """ Turbulence intensity profile at heights z for hub height zhub, 
        hub-height wind speed Vhub, and turbulence class turbc.
        
        Args:
            z (numpy array): array of heights in meters for velocity calcs
            zhub (float): hub-height of turbine in meters
            Vhub (float): hub-height mean velocity in m/s
            turbc (string): turbulence class from IEC 61400-1 Ed. 3
            
        Returns:
            Ti (numpy array): turbulence intensity at heights in z
    """

    
    Iref = IrefFromClass(turbc);                # reference hub-height TI  
    sigma1 = calcSigma1(Iref,Vhub)              # reference standard deviation
    V = VelProfile(z,zhub,Vhub);                # mean wind profile
    Ti = sigma1 / V;                            # TI profile
    
    return Ti


def IEC_Lambda1(zhub):
    """ Longitudinal scale parameter
    
        Args:
            zhub (float): hub-height of turbine in meters
        
        Returns:
            Lambda1 (float): longitudinal scale parameter
    """
    
    Lambda1 = (0.7 * zhub)*(zhub <= 60) \
        + 42*(zhub > 60);               # longitudinal scale parameter
    
    return Lambda1


def IEC_SpatialCoherence(zhub,Vhub,rsep,f):
    """ Spatial coherence function for hub-height wind speed Vhub and hub 
        height zhub at separation distance rsep and frequencies f.
        
        Args:
            zhub (float): hub-height of turbine in meters
            Vhub (float): hub-height mean velocity in m/s
            rsep (float): separation distance in m
            f (numpy array): frequencies in Hz
        
        Returns:
            Coh (numpy array): values of spatial coherence function
    """
    import numpy as np

    Lambda1 = IEC_Lambda1(zhub)         # longitudinal scale parameter
    Lc = 8.1*Lambda1;                   # coherence scale parameter

    Coh = np.exp(-12*np.sqrt( \
        np.power(f*rsep/Vhub,2) + \
        np.power(0.12*rsep/Lc,2) ));    # spatial coherence

    return Coh


def KaimalSpectrum(f,tau,sig):
    """ Kaimal spectrum (continuous, one-sided) for frequency f and time
        length scale tau = L/U.

        Args:
            f (numpy array): frequencies
            tau (float/numpy array): integral time scale (L/U)
            sig (float): standard deviation
       
        Returns:
            S (numpy array): Kaimal spectrum evaluated at f, tau, sig
    """
    import numpy as np

    S = (sig**2)*(4.*tau)/ \
        np.power(1.+6.*f*tau,5./3.);            # Kaimal 1972

    return S


def IEC_PSDs(zhub,Vhub,turbc,f):
    """ Continuous power spectral densities for all 3 turbulence components
    
        Args:
            zhub (float): hub-height of turbine in meters
            Vhub (float): hub-height mean velocity in m/s
            turbc (string): turbulence cass from IEC 61400-1 Ed. 3
            f (numpy array): frequencies in Hz
            
        Returns:
            Su (numpy array): longitudinal PSD
            Sv (numpy array): lateral PSD
            Sw (numpy array): vertical PSD
    """

    Iref = IEC_Iref(turbc);             # reference turbulence intensity
    sigma1 = IEC_Sigma1(Iref,Vhub);     # longitudinal standard deviation
    sigma2 = 0.8*sigma1;                # lateral standard deviaiton
    sigma3 = 0.5*sigma1;                # vertical standard deviation
    Lambda1 = IEC_Lambda1(zhub);        # longitudinal scale parameter
    L1 = 8.1*Lambda1;                   # longitudinal integral scale   
    L2 = 2.7*Lambda1;                   # lateral integral scale  
    L3 = 0.66*Lambda1;                  # vertical integral scale

    Su = KaimalSpectrum(f,L1/Vhub,sigma1);  # longitudinal spectrum
    Sv = KaimalSpectrum(f,L2/Vhub,sigma2);  # lateral spectrum
    Sw = KaimalSpectrum(f,L3/Vhub,sigma3);  # vertical spectrum

    return (Su,Sv,Sw)


def IEC_TC_simulation(turbc,zhub,Vhub,n_t,dt,rho,mu,y_grid,z_grid):
    """ Stochastic simulation of KSEC model from IEC with temporal
        coherence

        Args:

        Returns:
            u_grid (numpy array): [n_z x n_y x n_t] array of u(t)
    """
    import numpy as np
    
    # define variables to be subsequently used          
    n_y     = y_grid.size                       # no. of y-points in grid
    n_z     = z_grid.size                       # no. z-points in grid
    n_grid  = n_y*n_z                           # total no. grid points
    n_f     = uniqueComponents(n_t)          # no. unique frequencies
    df      = 1./(n_t*dt)                       # frequency resolution
    fs      = np.arange(n_f)*df                 # array of frequencies
    Y_grid, Z_grid = np.meshgrid(y_grid,z_grid) # 2D arrays of grid points
    y_vec   = Y_grid.reshape(n_grid)            # 1D array of grid points
    z_vec   = Z_grid.reshape(n_grid)            # 1D array of grid points

    # calculate matrix of distances for later spatial coherence
    DR = np.empty((n_grid,n_grid))
    for i in range(n_grid):
        for j in range(i,n_grid):
            dy  = y_vec[i] - y_vec[j]
            dz  = z_vec[i] - z_vec[j]
            DR[i,j] = np.sqrt(dy**2 + dz**2)
            DR[j,i] = np.sqrt(dy**2 + dz**2)
            
    # get array of magnitudes
    U_z = IEC_VelProfile(z_vec,zhub,Vhub)               # velocity profile
    Su = IEC_PSDs(zhub,Vhub,turbc,fs)[0]                # longitudinal PSD
    Xu_mags = np.empty((n_f,n_grid),dtype='complex')    # Fourier magintudes
    Xu_mags[:] = np.sqrt(Su*df/2).reshape(n_f,1)        # Fourier magnitudes
            
    # loop through frequencies, spatially correlating phases
    Xu_out      = np.empty((n_f,n_grid),dtype='complex')
    Xu_out[0,:] = U_z
    dphis_grid = wrappedCauchySample(n_f,n_grid,rho,mu)
    phis_grid  = wrap(np.cumsum(dphis_grid,axis=0))
    for i in range(1,n_f):
        
        # get Cholesky decomposition
        f = fs[i]
        Coh = IEC_SpatialCoherence(zhub,Vhub,DR,f)
        C = np.linalg.cholesky(Coh)
        
        # unspatially correlated Fourier component for grid
        phis = phis_grid[i,:].reshape(n_grid,1)
        Xu_noSC    =  Xu_mags[i,:].reshape(n_grid,1) * \
                        np.exp(1j*phis)
        
        # correlate
        Xu_SC      = np.dot(C,Xu_noSC)
        
        # save in array
        Xu_out[i,:] = np.squeeze(Xu_SC)

    # take IFT to get time series
    u_2D = np.fft.irfft(Xu_out,axis=0)*n_t
    u_grid = np.empty((n_z,n_y,n_t))
    i_tot = 0
    for i in range(n_z):
        for j in range(n_y):
            u_grid[i,j,:] = u_2D[:,i_tot]
            i_tot += 1

    return u_grid


def data2field(y_data,z_data,t_data,u_data,y_grid,z_grid,zhub=90.,Vhub=10.):
    """ Construct full turbulence field from set of discrete measurements
    
        Args:
            y_data (numpy array): y-locations of measurements
            z_data (numpy array): z-locations of measurements
            t_data (numpy array): measurement time stamps
            u_data (numpy array): (n_t x n_data) array of measurements
            y_grid (numpy array): y-positions of grid points
            z_grid (numpy array): z-positions of grid points
            
        Returns:
            u_grid (numpy array): (n_z x n_y x n_t) array of wind velocity
    """
    import sys
    libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
    if (libpath not in sys.path): sys.path.append(libpath)
    
    import JR_Library.main as jr
    import numpy as np
    from scipy import optimize
    
    # check array sizes
    if (y_data.size != z_data.size):
        raise ValueError('y_data and z_data must be the same size')
    if (y_data.size != u_data.shape[1]):
        raise ValueError('u_data must have the same number of columns as ' + \
            'the size of y_data')

    # define variables to be subsequently used          
    n_data  = y_data.size                       # no. of data points
    n_y     = y_grid.size                       # no. of y-points in grid
    n_z     = z_grid.size                       # no. z-points in grid
    n_grid  = n_y*n_z                           # total no. grid points
    n_all   = n_data + n_grid                   # all points in total
    n_t     = u_data.shape[0]                   # no. time steps
    dt      = (t_data[10]-t_data[0])/float(10)  # time step
    n_f     = jr.uniqueComponents(n_t)          # no. unique frequencies
    df      = 1./(n_t*dt)                       # frequency resolution
    fs      = np.arange(n_f)*df                 # array of frequencies
    Y_grid, Z_grid = np.meshgrid(y_grid,z_grid) # arrays of grid points
    y_all = np.concatenate((y_data,
                    Y_grid.reshape(n_grid)))    # vector of all y coordinates
    z_all = np.concatenate((z_data,
                    Z_grid.reshape(n_grid)))    # vector of all z coordinates
                    
    # calculate matrix of distances for later spatial coherence
    DR = np.empty((n_all,n_all))
    for i in range(n_all):
        for j in range(i,n_all):
            dy  = y_all[i] - y_all[j]
            dz  = z_all[i] - z_all[j]
            DR[i,j] = np.sqrt(dy**2 + dz**2)
            DR[j,i] = np.sqrt(dy**2 + dz**2)
            if ((i != j) and np.abs(DR[i,j] < 1e-10)):
                raise ValueError('Currently grid points and msmnt pts cannot be co-located.')
            
    # fit power law to measured wind
    U_z = np.mean(u_data,axis=0)
    fitfunc = lambda p, x: p[0] * x ** (p[1])
    errfunc = lambda p, x, y: (y - fitfunc(p, x))
    p_out = optimize.leastsq(errfunc, [10.,0.2],args=(z_data, U_z))[0]            
            
    # get spectral values of data
    U_data     = np.fft.rfft(u_data,axis=0)/n_t
    Umags_data = np.abs(U_data)
    
    # initialize array of spectral values, assign DC value from power law fit
    Umags_grid      = np.empty((n_f,n_grid))
    Umags_grid[0,:] = p_out[0]*np.power(z_all[n_data:],p_out[1])
    
    # set grid spectral values by finding/averaging values at closest height
    for i in range(n_grid):
        
        # find closest measurement height
        z_p = z_all[i]                          # grid point height
        dzs = np.abs(z_data-z_p)                # height diff
        i_closest = np.where(dzs == dzs.min())  # closest measurement index
        
        # average power of S(f)s at closest height to get grid S(f)
        S_closest = np.power( \
                    Umags_data[:,i_closest],2)   # S(f) of closest point(s)
        U_p = np.squeeze(np.sqrt( \
                    np.mean(S_closest,axis=1)))  # mean-power magnitude
        Umags_grid[:,i] = U_p                    # save magnitude
        
    # initialize output spectral matrix for ALL points (measurements + grid)
    U_all = np.empty((n_f,n_all),dtype='complex')
    
    # set DC values for data and grid points
    U_all[0,:n_data] = U_z
    U_all[0,n_data:] = p_out[0]*np.power(z_all[n_data:],p_out[1])
    
    # loop through frequencies, solving for phase pre-spatial correlation and
    #   then spatially correlating to get result
    for i in range(1,n_f):
        
        # get Cholesky decomposition
        f = fs[i]
        Coh = jr.IEC_SpatialCoherence(zhub,Vhub,DR,f)
        C = np.linalg.cholesky(Coh)
        
        # solve for un-spatially correlated data Fourier coefficients
        a = C[:n_data,:n_data]
        b = U_data[i,:].reshape(n_data,1)
        Xu_data = np.linalg.solve(a,b)
        
        # unspatially correlated Fourier component for grid
        phis_grid  = np.random.rand(n_grid,1) * 2 * np.pi
        Xu_grid    =  Umags_grid[i,:].reshape(n_grid,1) * \
                        np.exp(1j*phis_grid)
        
        # create entire vector
        Xu_all  = np.concatenate((Xu_data,Xu_grid),axis=0)
        
        # correlate
        Xs_all = np.dot(C,Xu_all)
        
        # save in array
        U_all[i,:] = np.squeeze(Xs_all)
        
        if (i == 1):
            Xu_grid_out = Xu_grid
    
    # take IFT to get time series
    u_all = np.fft.irfft(U_all,axis=0)*n_t
    u_grid = np.empty((n_z,n_y,n_t))
    i_tot = n_data
    for i in range(n_z):
        for j in range(n_y):
            u_grid[i,j,:] = u_all[:,i_tot]
            i_tot += 1
    u_data_out = u_all[:,:n_data]                   # verif: should = u_data
            
    out_dict = {}
    out_dict['u_grid']      = u_grid
    out_dict['u_data_out']  = u_data_out
    out_dict['p_out']       = p_out
    out_dict['Xu_grid_out'] = Xu_grid_out
            
    return out_dict
    

# ==============================================================================
# TURBSIM ANALYSIS
# ==============================================================================

def vectorSpatCoh(Xi,Xj,df):
    """ Takes [n_f] x [n_p] arrays Xi and Xj, where n_p is the number of pairs 
        of points to be used in the calculation and returns arrays of the 
        frequency and coherence.
        
        Args:
            Xi ([n_f x n_p] numpy array): Fourier coefficients at first point
            Xj ([n_f x n_p] numpy array): Fourier coefficients at second point
    """

    import numpy as np

    # generate frequency vector
    n_f = Xi.shape[0];
    f   = np.arange(0,n_f)*df;

    # calculate expected values
    EXiXi = np.real(np.mean(np.multiply(Xi,\
        np.conj(Xi)),axis=1));
    EXjXj = np.real(np.mean(np.multiply(Xj,\
        np.conj(Xj)),axis=1));
    EXiXj = np.mean(np.multiply(Xi,\
        np.conj(Xj)),axis=1);

    # calculate PSDs/CSDs
    Sii = 2 * EXiXi / df;
    Sjj = 2 * EXjXj / df;
    Sij = 2 * EXiXj / df;

    # calculate coherence
    Cohij = np.abs( np.divide( Sij, \
        np.multiply(np.sqrt(Sii),np.sqrt(Sjj)) ) );

    return (f, Cohij)


def TurbSimSpatCoh(fname,rsep):
    """ Loads turbsim file in fname and returns the spatial coherence at 
        frequency f given separation distance rsep.
        
        Args:
            fname (string): filename
            rsep (float): separation distance in meters
            
        Returns:
            
    """
    import pyts.io.main as io
    import numpy as np

    # read file
    tsout = io.readModel( fname );

    # extract grid
    grid = tsout.grid;

    # define useful parameters for simpler code
    n_p = grid.n_p;
    n_t = tsout.uhub.size;          # grid.n_t is not correct
    n_f = np.ceil((n_t-1)/2+1);     # unique components counting DC
    df  = grid.df;

    # initialize PSD/CSD arrays
    Xi_all = np.zeros( [n_f, 1]);
    Xj_all = np.zeros( [n_f, 1]);

    # loop through points in grid
    for ii in range(0,n_p):
        for jj in range(ii,n_p):

            # check if correct separation
            if ( grid.dist(ii,jj) == rsep ):

                # get grid indices
                [iiz,iiy] = grid.ind2sub(ii);
                [jjz,jjy] = grid.ind2sub(jj);

                # extract time histories
                ui = tsout.u[iiz,iiy,:];
                uj = tsout.u[jjz,jjy,:];

                # calculate Fourier coefficients
                Xi = np.fft.fft(ui)[0:n_f].reshape([n_f,1])/n_t;
                Xj = np.fft.fft(uj)[0:n_f].reshape([n_f,1])/n_t;

                # append to array
                Xi_all = np.hstack([Xi_all, Xi]);
                Xj_all = np.hstack([Xj_all, Xj]);

    # remove zeros from initial column
    Xi_all = Xi_all[:,1:];
    Xj_all = Xj_all[:,1:];

    # calculate spatial coherence of arrays
    f, Cohij = spatialCoherence(Xi_all,Xj_all,df);

    return (f, Cohij)


def TurbSimVelProfile(fname):
    """ An array of the mean velocities evaluated at each point in the TurbSim
        grid of file fname.
        
        Args:
            fname (string): TurbSim filename
            
        Returns:
            U (numpy array): mean velocities at heights in grid
    """
    import pyts.io.main as io
    import numpy as np

    # read file
    tsout = io.readModel(fname);

    # rename grid for simplicity
    grid = tsout.grid;

    # initialize U array
    U = np.empty( grid.shape );

    # loop through grid pts, evaluate U, and save
    for iy in range(0,grid.n_y):
	for iz in range(0,grid.n_z):
		U[iz,iy]  = np.mean( tsout.u[iz,iy,:] );

    return U


def TurbSimHHPSDs(fname):
    """ Frequency vector and hub-height spectra for TurbSim output.
    
        Args:
            fname (string): TurbSim filename
            
        Returns:
            Suk (numpy array): discrete longitudinal PSD
            Svk (numpy array): discrete lateral PSD
            Swk (numpy array): discrete vertical PSD
    """
    import pyts.io.main as io

    # read file
    tsout = io.readModel(fname);

    # calculate PSDs
    Suk = signal2Sk(tsout.uhub)
    Svk = signal2Sk(tsout.vhub)
    Swk = signal2Sk(tsout.whub)
    
    return (Suk,Svk,Swk)

    
def TurbSimHHPDDs(fname):
    """ Hub-height phase differences
    """
    import pyts.io.main as io
    import numpy as np

    # read file
    tsout = io.readModel(fname);

    n_t = tsout.uhub.size;          # stored grid.n_t incorrect
    n_f = np.ceil((n_t-1)/2+1);     # unique Fourier (counts DC)

    # Fourier transform
    Xu = np.fft.fft( tsout.uhub )[0:n_f]/n_t;
    Xv = np.fft.fft( tsout.vhub )[0:n_f]/n_t;
    Xw = np.fft.fft( tsout.whub )[0:n_f]/n_t;

    # Phase differences
    dthetau = np.angle( np.divide(Xu[1:], Xu[:-1]) );
    dthetav = np.angle( np.divide(Xv[1:], Xv[:-1]) );
    dthetaw = np.angle( np.divide(Xw[1:], Xw[:-1]) );

    return (dthetau,dthetav,dthetaw)


# ==============================================================================
# MAPPINGS
# ==============================================================================

def X2Sk(X):
    """ Convert Fourier array to discrete, one-sided spectral mass values.
        
        Args:
            X (numpy array): Fourier vector
        
        Returns:
            Sk (numpy array): array of discrete spectral values
    """
    import numpy as np

    n_t = X.size;                           # no. elements
    n_f = uniqueComponents(n_t);            # no. unique comps
    
    Xuniq = X[:n_f];                        # unique Fourier comps
    
    Sk = np.real(2*np.multiply(Xuniq, \
            np.conj(Xuniq)));               # discrete spectrum

    return Sk


def Sk2Xuniq(Sk):
    """ Convert array of discrete spectral values in Sk to array of unique 
        Fourier magnitudes X.

        Args:
            Sk (numpy array): array of discrete spectral values
        
        Returns:
            Xuniq (numpy array): array of unique Fourier components
    """
    import numpy as np

    Xuniq = np.sqrt(Sk/2);

    return Xuniq


def Xuniq2X(Xuniq,n_t):
    """ Convert array of unique Fourier components to full Fourier array.
        
        Args:
            Xuniq (numpy array): array of unique Fourier components
            n_t (int): no. of components of full array
        
        Returns:
            X (numpy array): array of all Fourier components
    """
    import numpy as np

    if (n_t % 2):                           # even no. of elements
        Xuniq[0] = np.abs(Xuniq[0])         # ensure real signal
        X = np.concatenate((Xuniq,np.conj( \
            Xuniq[1:][::-1])),axis=0);
    else:                                   # odd no. of elements
        Xuniq[0] = np.abs(Xuniq[0])         # ensure real signal
        Xuniq[-1] = np.abs(Xuniq[-1])       # ensure real signal
        X = np.concatenate((Xuniq,np.conj( \
                    Xuniq[1:-1][::-1])),axis=0);

    return X


def signal2Sk(x):
    """ Return array of discrete spectrum for signal x.

        Args:
            x (numpy array): array of time history
       
        Returns:
            Sk (numpy array): array of discrete spectral values
    """
    import numpy as np

    n_t  = x.shape[0];                      # no. of elements
    X    = np.fft.fft(x,axis=0)/n_t;        # Fourier vector
    Sk  = X2Sk(X);                          # discrete PSD

    return Sk


def samplePhaseCoherence(theta,axis=-1):
    """ Return concentration and location parameters for a sample of wrapping
        random variables.
        
        Args:
            theta (numpy array): 1D array of sample of angles
       
        Returns:
            rho (float): concentration parameter
            mu (float): location parameter
    """
    import numpy as np
        
    z = np.exp(1j*theta)
    V = np.mean(z,axis=axis)
    rho = np.abs(V)
    mu  = np.angle(V)
        
    return (rho,mu)


def signalPhaseDifferences(x,axis=-1):
    """ Return phase differences for a time history.
    
        Args:
            x (numpy array): 1D array of time history
       
        Returns:
            rho (float): concentration parameter
            mu (float): location parameter
    """
    import numpy as np
    
    # if (2+)D array is fed in, halt with error
#    if ((len(x.shape)>1) and (x.shape[0] != 1 and x.shape[1] != 1)):
#        errStr = 'signalPhaseDifferences only works on 1D arrays'
#        raise AttributeError(errStr)
    
    Xuniq = np.fft.rfft(x,axis=axis)            # Fourier vector
    thetas = np.angle(Xuniq)                    # Fourier angles
    dtheta = wrap(np.diff(thetas,axis=axis))    # phase differences
    
    return dtheta


def signalPhaseCoherence(x):
    """ Return concentration and location parameters for a time history.
    
        Args:
            x (numpy array): 1D array of time history
       
        Returns:
            rho (float): concentration parameter
            mu (float): location parameter
    """
    import numpy as np
    
    # if (2+)D array is fed in, halt with error
    if ((len(x.shape)>1) and (x.shape[0] != 1 and x.shape[1] != 1)):
        errStr = 'signalPhaseCoherence only works on 1D arrays'
        raise AttributeError(errStr)

    # bypass if all nans
    if np.all(np.isnan(x)):
        rho, mu = np.NAN, np.NAN

    # otherwise remove nans and proceed
    else:
        
        nan_idx = np.isnan(x)                       # nan indices
        x_nonan = x[np.logical_not(nan_idx)]        # remove nanes
        dtheta  = signalPhaseDifferences(x_nonan)   # phase differences
        rho, mu = samplePhaseCoherence( dtheta)     # temp coh parameters
    
    return (rho,mu)


def NRELfname2time(fname):
    """ Convert filename to float representing time stamp

        Args:
            fname (string): name of file

        Returns:
            timestamp (float): float representing timestamp
    """

    # extract date information
    year     = int(fname[6:10])
    month    = int(fname[0:2])
    day      = int(fname[3:5])
    hour     = int(fname[11:13])
    minute   = int(fname[14:16])
    time_tup = (year,month,day,hour,minute)

    # convert tuple to float
    time_flt = timetup2flt(time_tup)

    return time_flt


def time2fpath(dataset,timestamp):
    """ Convert float to path to data structure

        Args:
            timestamp (float or tuple): float or tuple (year,month,day,hour,minute)
                                        representing record start time

        Returns:
            fpath (string): path to data 
            
    """
    import time
    import os
    import glob

    if (dataset in ['NREL','fluela']):

        # get tuple of timestampe
        if isinstance(timestamp,tuple):
                time_tup = timestamp
        elif isinstance(timestamp,(int,float)):
                time_tup = timeflt2tup(timestamp)
        else:
            errStr = 'Invalid {} for timestamp'. \
                         format(type(timestamp))
            raise TypeError(errStr)

        # convert tuple to string values
        yearS  = str(time_tup[0])
        monthS = str(time_tup[1]).zfill(2)
        dayS   = str(time_tup[2]).zfill(2)
        hourS  = str(time_tup[3]).zfill(2)
        minS   = str(time_tup[4]).zfill(2)

        # get directory path
        basedir = getBasedir(dataset)
        dirpath = os.path.join(basedir,yearS,monthS,dayS)
        
        # get filename
        fname_part = '_'.join([monthS,dayS,yearS,hourS,minS])
        fpath = glob.glob(os.path.join(dirpath,fname_part)+'*.mat')

        # check if file doesn't exist
        if len(fpath) > 0:
            fpath = fpath[0]
        else:
            errStr = 'File {} does not exist.'.format(fname_part)
            raise IOError(errStr)

    else:
        errStr = 'Dataset {} is not coded.'.format(dataset)
        raise ValueError(errStr)
    
    return fpath


def timeflt2tup(time_flt):
    """ Convert time in float to tuple

        Args:
            time_flt (float): timestamp in float format

        Returns:
            time_tup (tuple): (year,month,day,hour,minute)
    """
    import time

    time_tup = time.gmtime(time_flt)[:5]    # convert float to tuple

    return time_tup


def timetup2flt(time_tup):
    """ Convert time in float to tuple

        Args:
            time_tup (tuple): (year,month,day,hour,minute)
            
        Returns:
            time_flt (float): timestamp in float format
            
    """
    import calendar

    time_tup = time_tup[:5]                 # ignore later entries in tuple
    time_tup += (0,0,0,)                    # seconds, milliseconds, zre zero
    time_flt = calendar.timegm(time_tup)    # convert to float

    return time_flt


def timeUTC2local(dataset,UTCtime_flt):
    """ Convert timestamp in UTC (float) to timestamp
        in local time for dataset (float)

        Args:
            dataset (string): flag to indicate dataset
            UTCtime_flt (int/flt/numpy array): time or array of times
                in UTC

        Returns:
            loctime_flt (int/flt/numpy array): time or array of times
                in local time
    """
    import datetime
    import numpy as np

    # convert integers/floats to 1D numpy arrays for calculations
    if (isinstance(UTCtime_flt,(int,float))):
        UTCtime_flt = np.array([UTCtime_flt])

    # initialize local time output
    loctime_flt = np.empty(UTCtime_flt.shape)
    
    if dataset == 'NREL':

        # convert each timestamp in a loop
        for i in range(loctime_flt.size):

            UTCtime_tup = timeflt2tup(\
                UTCtime_flt[i])                 # UTC tuple
            UTCtime_dat = datetime.datetime(\
                UTCtime_tup[0],UTCtime_tup[1],\
                UTCtime_tup[2],UTCtime_tup[3],\
                UTCtime_tup[4])                 # UTC datetime object
            time_delta  = datetime.timedelta(\
                hours = -7)                     # time difference
            loctime_dat = UTCtime_dat \
                          + time_delta          # local datetime object
            loctime_flt[i] = timetup2flt(\
                loctime_dat.timetuple()[:5])    # local float

    else:
        errStr = 'Dataset \"{}\" not coded.'.format(dataset)
        raise KeyError(errStr)

    return np.squeeze(loctime_flt)


def timeflt2arr(time_flt):
    """ Expand timestamp in float to [year,month,day,hour,minute]

        Args:
            time_flt (int/flt/numpy array): N timestamps in float form

        Returns:
            time_vec (numpy array): Nx5 array of timestamps in
                "expanded" format -- [year,month,day,hour,minute]
    """
    import numpy as np

    # convert integers/floats to 1D numpy arrays for calculations
    if (isinstance(time_flt,(int,float))):
        time_flt = np.array([time_flt])

    # initializze output array
    time_vec = np.empty((time_flt.size,5))

    # iteratively convert each timestamp
    for i in range(time_flt.size):
        time_tup = timeflt2tup(time_flt[i])
        time_vec[i,:] = np.asarray(time_tup)

    return time_vec


def field2datfield(dataset,field,ht):
    """ Dataset-specific fieldname corresponding with custom fieldname
        and height

        Args:
            dataset (string): flag for dataset
            field (string): custom fieldname
            ht (float): measurement height

        Returns:
            datfield (string): dataset-specific fieldname
    """

    ht = int(ht)                        # height must be int for string calculations

    if dataset == 'NREL':

        if   (field == 'Sonic_u'):        datfield = 'Sonic_u_' + str(ht) + 'm'
        elif (field == 'Sonic_v'):        datfield = 'Sonic_v_' + str(ht) + 'm'
        elif (field == 'Sonic_w'):        datfield = 'Sonic_w_' + str(ht) + 'm'
        elif (field == 'Sonic_T'):        datfield = 'Sonic_Temp_rotated_' + str(ht) + 'm'
        elif (field == 'Temperature'):    datfield = 'Air_Temp_' + str(ht) + 'm'
        elif (field == 'Wind_Speed_Cup'): datfield = 'Cup_WS_' + str(ht) + 'm'
        elif (field == 'Wind_Direction'): datfield = 'Vane_WD_' + str(ht) + 'm'
        elif (field == 'Precipitation'):  datfield = 'PRECIP_INTEN'
        elif (field == 'Dewpt_Temp'):     datfield = 'Dewpt_Temp_' + str(ht) + 'm'
        elif (field == 'Pressure'):       datfield = 'Baro_Presr_3m'
        else:
            errStr = 'Unknown custom field {} for dataset \"{}\"'.format(field,dataset)
            raise AttributeError(errStr)
            
        # check that datafield is valid
        if (not check_datfields(dataset,datfield)):
            raise ValueError('Invalid height {} for field {}'.format(\
                ht,field))

    elif dataset == 'fluela':

        if   (field == 'Sonic_u'):         datfield = 'Sonic_u_' + str(ht) + 'm'
        elif (field == 'Sonic_v'):         datfield = 'Sonic_v_' + str(ht) + 'm'
        elif (field == 'Sonic_w'):         datfield = 'Sonic_w_' + str(ht) + 'm'
        elif (field == 'Sonic_T'):         datfield = 'Sonic_Temp_rotated_' + str(ht) + 'm'
        elif (field == 'Sonic_Cup'):       datfield = 'Sonic_CupEqHorizSpeed_' \
                                                         + str(ht) + 'm'
        elif (field == 'Sonic_Direction'): datfield = 'Sonic_direction_' \
                                                         + str(ht) + 'm'
        else:
            errStr = 'Unknown custom field {} for dataset \"{}\"'.format(field,dataset)
            raise AttributeError(errStr)
            
        # check that datafield is valid
        if (not check_datfields(dataset,datfield)):
            raise ValueError('Invalid height {} for field {}'\
                             .format(ht,field))

    else:
        errStr = 'Dataset \"{}\" not coded.'.format(dataset)
        raise KeyError(errStr)

    return datfield


def grid2ticklist(grid_locs):
    """ Convert list of grid locations to tuple array of tick positions and
        and labels for plotting TurbSim plots

        Args:
            grid_locs (numpy ndarray): grid locations

        Returns:
            ticks (list): list of tuples (tick_loc,tick_label)
    """
    n_grid = grid_locs.size                     # number grid points
    i_lab  = [0,n_grid/2,n_grid-1]              # indices w/text labels

    ticks  = []
    for i in range(n_grid):
        if (i in i_lab):                        # add labels 1st, last, mid
            ticks.append((grid_locs[i],'{:.0f}'.format(grid_locs[i])))
        else:                                   # no labels elsewhere
            ticks.append((grid_locs[i],''))

    return ticks



def C2K(T_c):
    """ Celsius to Kelvin
    
        Args:
            T_c (float,numpy array): temperature in Celsius
            
        Returns:
            T_k (float, numpy array): temperature in Kelvin
    """
    return T_c + 273.15
    

def K2C(T_k):
    """ Kelvin to Celsius
    
        Args:
            T_k (float, numpy array): temperature in Kelvin
            
        Returns:
            T_c (float,numpy array): temperature in Celsius
    """
    return T_k - 273.15


# ==============================================================================
# MISCELLANEOUS
# ==============================================================================

def wrap(theta):
    """ Wrap angle to [-pi,pi).
    
        Args:
            theta (float): angle in radians
    
        Returns:
            theta_out (float): wrapped angle
    """
    import numpy as np

    theta_out = theta % (2*np.pi)           # mod to [0,2pi)
    theta_out = np.multiply(theta_out<=np.pi,theta_out) + \
        np.multiply(theta_out>np.pi, \
                    theta_out-2*np.pi)      # subtract 2pi from [pi,2pi)
    
    return theta_out


def uniqueComponents(n_t):
    """ Number of unique Fourier components for a vector of length n_t.

        Args:
           n_t (int): vector length
        
        Returns:
           n_f (int): no. unique components
    """
    import numpy as np

    n_f = int(np.ceil((float(n_t)-1) \
        /2)+1);

    return n_f


def spectralScale(Sk,sig,n_t):
    """ Spectral scaling factor for discrete power spectral coefficients Sk.
        
        Args:
            Sk (numpy array): array of discrete spectral values
            sig (float): desired standard deviation of time history
            n_t (int): number of components in time history
        
        Returns:
            alpha (float): spectral scaling factor
    """
    import numpy as np
    
    # if (2+)D array is fed in, halt with error
    if ((len(Sk.shape)>1) and (Sk.shape[0] != 1 and Sk.shape[1] != 1)):
        print 'ERROR: spectralScale only works on 1D arrays'
        return []

    Xuniq   = Sk2Xuniq(Sk)                      # unique Fourier components
    Xuniq[0] = 0                                # neglect DC component
    X       = Xuniq2X(Xuniq,n_t)                # entire Fourier vector
    sumXkSq = np.sum(np.power(np.abs(X),2));    # sum |X[k]|^2
    alpha   = np.sqrt((n_t-1.)/n_t \
        *(sig**2)/sumXkSq);                     # scaling factor
    
    return alpha


def wrappedCauchyPDF(theta,rho,mu):
    """ Probability density function of wrapped Cauchy distribution.
    
        Args:
            theta (numpy array): angles (rads) to evaluate PDF
            rho (float): concentration parameter
            mu (float): location parameter
       
        Returns:
            f (numpy array): evaluation of wrapped Cauchy PDF at angles in theta
    """
    import numpy as np
    
    if ((rho < 0) or (rho>1)):
        print 'rho must be between 0 and 1'
        return []
    
    f = 1/(2*np.pi)*(1 - rho**2) / \
        (1 + (rho**2) - 2*rho*np.cos(theta-mu))
    
    return f


def removeSpines(ax):
    """ Remove the top and right spines on a Python plot

        Args:
            ax (axis handle): handle to axes to remove spines
    """
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    return


def nandetrend(x,y):
    """ Linear detrend ignoring nan values in y (removes bias too)

        Args:
            x (numpy array): x values
            y (numpy array): y values

        Returns:
            y_det (numpy array): detrended values with nans unchanged
    """
    import numpy as np
    
    # check that x and y are the same size
    if (x.size != y.size):
        errStr = 'Sizes for x and y are not equal ({} and {})' \
            .format(x.size,y.size)
        raise ValueError(errStr)

    # remove singleton dimensions if any
    if (len(x.shape) > 1):
        x = np.squeeze(x)
    if (len(y.shape) > 1):
        y = np.squeeze(y)

    # initialize output
    y_det   = np.empty(y.shape)

    # check if all NaNs
    if np.all(np.isnan(y)):
        y_det[:] = np.NAN

    # linear detrend if 1+ non-nan values
    else:

        # find NaN indices of y
        nan_idx = np.isnan(y)               

        # isolate non-NaN values
        xnew = x[np.logical_not(nan_idx)]                        
        ynew = y[np.logical_not(nan_idx)]

        # determine least-squares linear regression
        A = np.vstack([xnew, \
                       np.ones(len(xnew))]).T
        m, c = np.linalg.lstsq(A, ynew)[0]

        # remove trend
        y_det = y - m*x - c
    
    return y_det


def remove_spikes(x, spikeWidth=6, P=0.9, beta=10.):
    """ Linearly interpolate spikes in time series

        Spikes are detected by 1) finding standard deviation of the
        bottom P% differences, 2) finding ``extremal'' values (above
        beta times the standard deviation of the bottom P% differences),
        3) finding extremal values withing spikeWidth points of one
        another.

        Args:
            x (numpy array): time series
            spikeWidth (int): optional, max width of spike (no. of points)
            P (float): optional, threshold percentage for differences
            beta (float): optional, scaling factor for spike detection

        Returns:
            x_cl (numpy array): time series with spikes removed
    """
    import numpy as np

    N        = x.size                           # length of record
    x_cl     = np.copy(x)                       # cleaned version
    n_spikes = 0                                # initialize no. of spikes
    
    diffs    = x[1:] - x[:-1]                   # velocity differences
    absdiffs = abs(diffs)                       # sorted diff magnitudes
    thresh   = np.sort(absdiffs)[int(N*P)]      # absdiff P quantile value
    diffs_P  = diffs[np.where(abs(diffs) <= \
                    thresh)]                    # diffs below P% threshold
    sig_x    = np.std(diffs_P)                  # std dev of P% differences
    ext_idx  = np.where(absdiffs > \
        beta*sig_x)[0]                          # extremal indices

    # remove spike at end of record
    if np.any(ext_idx >= N-spikeWidth):
        spike_idx = ext_idx[np.where( \
            ext_idx >= N-spikeWidth)[0][0]]     # get spike index
        x_cl[spike_idx+1:] = x_cl[spike_idx]    # append clean value
        ext_idx = ext_idx[np.where( \
            ext_idx < N-spikeWidth)]            # update extremal indices
        n_spikes += 1
            
    # remove spike at beginning of record
    if np.any(ext_idx <= spikeWidth):
        spike_idx = ext_idx[np.where( \
            ext_idx <= spikeWidth)[0][-1]]      # get spike index
        x_cl[:spike_idx+1] = x_cl[spike_idx+1]  # append clean value
        ext_idx = ext_idx[np.where( \
            ext_idx > spikeWidth)]              # update extremal indices
        n_spikes += 1
            
    # loop through and remove spikes in middle of record
    i = 0
    while (i < ext_idx.size-1):
        
        # calculate distance between extremal values
        idx_left, idx_right = ext_idx[i], ext_idx[i+1]
        ext_dst = idx_right - idx_left
        
        # detect spike if two extremal values with opp sign are near each other
        spike = (ext_dst < spikeWidth) and \
            (diffs[idx_left]*diffs[idx_right] < 0)
            
        # linearly interpolate if spike is detected
        if spike:
            xp = np.array([idx_left,idx_right+1])
            yp = np.array([x_cl[idx_left],x_cl[idx_right+1]])
            x_cl[idx_left+1:idx_right+1] = \
                np.interp(np.arange(idx_left+1,idx_right+1),xp,yp)
            n_spikes += 1
        
        # continue looping through extremal values
        i += 1        
        
    return x_cl, n_spikes


def cleantimeseries(t,x,spikeWidth=6, P=0.9, beta=10.):
    """ Remove spikes and remove linear trend (same mean)

        Args:
            t (numpy array): time steps
            x (numpy array): raw time series

        Returns:
            x_cl (numpy array): time series with spikes
                            removed and detrended
    """
    import numpy as np

    # remove spikes
    x_nospike, n_spikes = remove_spikes(x,spikeWidth,P,beta)

    # remove linear trend
    x_cl = nandetrend(t,x_nospike) + np.nanmean(x_nospike)

    return x_cl, n_spikes


def is_quantized(x):
    """ Check if time series x is quantized
    
        Args:
            x (numpy array): time series
            
        Returns:
            flag (boolean): flag if quantized
    """
    import numpy as np
        
    # parameters to detect quantization
    n_zero    = 4                   # number of allowable consecutive zeros
    quant_tol = 1e-12               # tolerance for defining zero
    n_steps   = 10                  # no. quantization steps before flagging

    diffs      = x[1:] - x[:-1]     # point-to-point differences
    zero_jumps = diffs < quant_tol  # differences that are zero

    quant = np.zeros(n_zero)        # subarray for searching for quants

    # loop through elements of x searching for quants
    quant_idx = []
    for i in range(diffs.size-n_zero+1):
        if np.all(diffs[i:i+n_zero] == quant):
            quant_idx.append(i)

    if (len(quant_idx) > n_steps):
        flag = 1
    else:
        flag = 0

    return flag
    
def sheared_axes(fig,rect,x_ticks,y_ticks,skew_ang=3.14159/8):
    """ Define floating axes and sheared axes for plotting sheared data
    
        Args:
            fig (matplotlib figure): handle to matplotlib figure
            rect (int): 3-digit integer indicating subplots (e.g., 121)
            x_ticks (list): list of tuples (x_tick,x_ticklabel)
            y_ticks (list): list of tuples (y_tick,y_ticklabel)
            skew_ang (float,opt): skew angle in radians
            
        Returns:
            ax_flt (floating subplot): handle to floating axes
            ax_shr (parasite subplot): handle to sheared axes
    """
    from matplotlib.transforms import Affine2D
    import mpl_toolkits.axisartist.floating_axes as floating_axes
    from mpl_toolkits.axisartist.grid_finder import FixedLocator, DictFormatter
    
    # define affine transformation
    tr = Affine2D().skew(0,skew_ang)
    
    # get locations and labels of xticks
    xtick_locations = FixedLocator([xtick[0] for xtick in x_ticks])
    xtick_labels    = DictFormatter(dict(x_ticks))

    # get locations and labels of yticks
    ytick_locations = FixedLocator([ytick[0] for ytick in y_ticks])
    ytick_labels = DictFormatter(dict(y_ticks))
    
    # get extremes from tick labels
    x_lo, x_hi = x_ticks[0][0], x_ticks[-1][0]
    y_lo, y_hi = y_ticks[0][0], y_ticks[-1][0]
    
    # define grid helper to transform between data coors and affine coors
    grid_helper = floating_axes.GridHelperCurveLinear(tr,
                                extremes=(x_lo, x_hi, y_lo, y_hi),
                                grid_locator1=xtick_locations,
                                grid_locator2=ytick_locations,
                                tick_formatter1=xtick_labels,
                                tick_formatter2=ytick_labels,
                                )

    # create the floating axes
    ax_flt = floating_axes.FloatingSubplot(fig, rect, \
        grid_helper=grid_helper)
        
    # add floating axis to figure
    fig.add_subplot(ax_flt)

    # define rotated axes
    ax_shr = ax_flt.get_aux_axes(tr)

    return ax_flt, ax_shr


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    """
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.

    Example
    -------
        >> orig_cmap = matplotlib.cm.RdBu_r
        >> new_cmap  = jr.shiftedColorMap(orig_cmap)
          
    Found: http://stackoverflow.com/questions/7404116/...
            defining-the-midpoint-of-a-colormap-in-matplotlib
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = mpl.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap
