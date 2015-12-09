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
    - FAST Analysis:
        Analyze FAST results, create turbine models
    - Mappings:
        Time UTC to local, general field to dataset-specific, etc.
    - Miscellaneous:
        Other functions

Restructured 2015-06-30
Jenni Rinker, Duke University
"""

# used modules
import numpy as np
import scipy.io as scio

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

def loadtimeseries(dataset,field,ID,data_in):
    """ Load the time series for a given dataset, field and timestamp 
        or structure

        Args:
            dataset (string): flag to indicate dataset
            field (string): general field name
            ID (int): device identifier (often mstmnt height)
            data_in (float/tuple or dictionary): float or tuple representing 
                time value or high-frequency data structure

        Returns:
            outdict (dictionary): keys = ['raw','clean','flags']
    """
        
    # make sure datfield is valid before proceeding
    datfield = field2datfield(dataset,field,ID)
    if (not check_datfields(dataset,datfield)):
        raise KeyError('Invalid datafield \"{}\".'.format(\
                    datfield))

    # turn off warnings about invalid entries
    np.seterr(invalid='ignore')
    
    # define max number of spikes
    n_spikes_max = 20

    if (dataset in ['NREL','fluela']):

        # set data ranges, desired length of time series
        dataRng = dataRanges(dataset,datfield)
        N, dt   = datasetSpecs(dataset)[:2]
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
        
    elif (dataset == 'CM06'):

        # set data ranges, desired length of time series
        N, dt   = datasetSpecs(dataset)[:2]
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
            x_raw = np.squeeze(struc[datfield])            
            # flag 5003: data are all NaNs
            if (np.all(np.isnan(x_raw))):
                flags.append(5003)                
        except:
            # flag 1006: data are not in high-frequency structure
            x_raw = np.empty(N)
            x_raw[:] = np.nan
            flags.append(1006)
            
        # copy raw data for analysis, make sure is float for NaN values
        x_cl  = x_raw.astype(float)
                
        # process time series if not group warning
        if ('grp' not in datfield):     
            
            # get range of acceptable data
            dataRng = dataRanges(dataset,datfield)           
                
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
    
            # remove spikes/detrend sonics, add DC gain
            # detrend everything else
            if (len(flags) == 0):
                if any([key in datfield for key in ('Ux','Uy','Uz','Ts')]):
                    # remove spikes, linear detrend
                    x_cl, n_spikes = cleantimeseries(t,x_cl)
                    
                    # subtract dc_gain
                    x_cl = x_cl - dcgain(dataset,field,ID)
                    
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


def mygenfromtxt(fname,header=True,units=False,delimiter='\t'):
    """ Load delimited array into numpy with optional header and units. Assumes
        units line will only be present if header line present.
    
        Args:
            fname (string): filename or path to file
            header (Boolean): whether to read header line
            units (Boolean): whether to read units line
            delimiter (string): string used to separate values
            
        Returns:
            data (numpy array): numeric data read from the text file
            header (list): strings from header line (if header=True)
            units (list): strings from units line (if header,units=True)
    """
    
    
    # read header and unit lines if appropriate
    with open(fname,'r') as f:
        
        # read header line if present
        if header:
            header_text = f.readline().strip('\n').split('\t')
            
            # read units line if present
            if units:
                units_text = f.readline().strip('\n').split('\t')
            else:
                units_text = []
                
        else:
            header_text = []
            
    # use numpy to load numeric data
    data = np.genfromtxt(fname,delimiter=delimiter,
                         skip_header=header+units)
    
    # return different values based on inputs
    if not header:
        return data
    elif (header and (not units)):
        return data, header_text
    else:
        return data, header_text, units_text


def ReadFASTFile(fname):
    """ Read FAST data into numpy array (v7.02)
    """
    
    
    n_skip = 8
    FASTDict = {}
    
    with open(fname,'r') as f:
        for i in range(6):
            f.readline()
        Fields = f.readline().strip('\n').split()
        Units  = f.readline().strip('\n').split()
    Data = np.genfromtxt(fname,skip_header=n_skip,delimiter='\t')
    
    # save everything in the dictionary
    FASTDict['Data']   = Data
    FASTDict['Fields'] = Fields
    FASTDict['Units']  = Units
    
    return FASTDict


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
              'MO_Length_interp','MO_Length_near','MO_Length_virt']
    elif (dataset == 'fluela'):
        fields = ['Record_Time','Processed_Time','Height','Sonic_Cup', \
              'Sonic_Direction', 'Mean_Wind_Speed', \
              'Sigma_u','Concentration_u','Location_u', \
              'Sigma_v','Concentration_v','Location_v', \
              'Sigma_w','Concentration_w','Location_w', \
              'up_wp','vp_wp','wp_Tp','up_vp','tau_u','tau_v','tau_w', \
              'MO_Length','MO_Length_virt']
    elif (dataset == 'CM06'):
        fields = ['Record_Time','Processed_Time','ID', \
              'Sonic_Cup', 'Sonic_Direction', 'Specific_Humidity',\
              'Mean_Wind_Speed', \
              'Sigma_u','Concentration_u','Location_u', \
              'Sigma_v','Concentration_v','Location_v', \
              'Sigma_w','Concentration_w','Location_w', \
              'up_wp','vp_wp','wp_Tp','up_vp','Tau_u','Tau_v','Tau_w', \
              'MO_Length','MO_Length_virt','Tbar_K','Tbar_K']
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)
        
    return fields


def check_datfields(dataset,datfield):
    """ Define list of fields of all possible correct raw fieldnames for 
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
    elif (dataset == 'CM06'):
        fields = ['Uy_1_3(1,1)', 'Ux_1_1(1,1)', 'Uz_2_3(1,1)', 'Ux_3_3(1,1)', \
        'rec_numbr_mas_1', 'Ux_2_2(1,1)', 'Ux_3_4(1,1)', 'Uz_2_4(1,1)', \
        'Ts_3_4(1,1)', 'Uy_2_3(1,1)', 'Ux_1_3(1,1)', 'Uy_2_4(1,1)', \
        'Ts_2_1(1,1)', 'Ts_2_2(1,1)', 'Uy_1_4(1,1)', 'Uz_1_2(1,1)', \
        'grp_warning_1(1)', 'Uz_3_2(1,1)', 'Uz_3_1(1,1)', 'Ux_2_3(1,1)', \
        'Uy_1_2(1,1)', 'Uy_3_3(1,1)', 'Ts_1_2(1,1)', 'rec_numbr_slv_3', \
        'rec_numbr_slv_2', 'Uy_2_2(1,1)', 'Ts_2_4(1,1)', 'Ts_3_1(1,1)', \
        'Ux_1_4(1,1)', 'Uz_1_3(1,1)', 'Uz_3_4(1,1)', 'Ux_3_1(1,1)', \
        'Ts_1_4(1,1)', 'Uz_2_1(1,1)', 'Uy_3_4(1,1)', 'Uz_1_1(1,1)', \
        'kh2o_1_1', 'Uy_2_1(1,1)', 'grp_warning_3(1)', 'Ux_1_2(1,1)', \
        'Ts_1_1(1,1)', 'Uy_1_1(1,1)', 'Uy_3_1(1,1)', 'grp_warning_2(1)', \
        'Ts_2_3(1,1)', 'TIMESTAMP', 'Uy_3_2(1,1)', 'Ts_3_2(1,1)', \
        'Ux_2_4(1,1)', 'kh2o_2_1', 'Ux_2_1(1,1)', 'Ts_1_3(1,1)', 'Ts_3_3(1,1)', \
        'Ux_3_2(1,1)', 'Uz_3_3(1,1)', 'Uz_1_4(1,1)', 'Uz_2_2(1,1)']

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
            IDs (numpy array): anemometer IDs (height or number)
    """
    

    if (dataset == 'NREL'):
        n_t     = 12000
        dt      = 0.05
        IDs     = np.array([15,30,50,76,100,131])
    elif (dataset == 'fluela'):
        n_t     = 6000
        dt      = 0.10
        IDs     = np.array([36,54,75])
    elif (dataset == 'CM06'):
        n_t     = 12000
        dt      = 0.05
        IDs     = range(1,13)
    elif (dataset == 'texastech'):
        n_t     = 30000
        dt      = 0.02
        IDs     = np.array([1,2,4,10,17,47,75,116,158,200])
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return n_t, dt, IDs


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
              
    elif (dataset == 'CM06'):

        # define data ranges
        #    Sonic data ranges taken from CSAT 3 specifications
        #    Hygrometer data from KH20 specifications
        if ('Ux' in datfield):
            dataRng = [-29.95,29.95]
        elif ('Uy' in datfield):
            dataRng = [-29.95,29.95]
        elif ('Uz' in datfield):
            dataRng = [-7.95,7.95]
        elif ('Ts' in datfield):
            dataRng = [-29.95,49.95]
        elif ('kh2o' in datfield):
            dataRng = [  1.75,19.25] 
        else:
              raise KeyError('Field {} not recognized.'.format(datfield))

    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return dataRng


def getBasedir(dataset,drive='G:'):
    """ Get path to base directory and check if it exists

        Args:
            dataset (string): flag for dataset
            drive (string): drive identifier [opt]

        Returns:
            basedir (str): path to top level of data directory
    """
    import os
    import platform

    if (dataset == 'NREL'):
        if (platform.system() == 'Linux'):
            basedir = '/media/jrinker/JRinker SeaGate External/' + \
                        'data/nrel-20Hz/'
        elif (platform.system() == 'Windows'):
            basedir = drive + '\\data\\nrel-20Hz'
        if not os.path.exists(basedir):
            errStr = 'Incorrect or unavailable base ' + \
                     'directory for dataset \"{}\".'.format(dataset)
            raise IOError(errStr)
    elif (dataset == 'fluela'):
        if (platform.system() == 'Linux'):
            basedir = '/media/jrinker/JRinker SeaGate External/' + \
                        'data/fluela-high_freq/'
        elif (platform.system() == 'Windows'):
            basedir = drive + '\\data\\fluela-high_freq'
        if not os.path.exists(basedir):
            errStr = 'Incorrect or unavailable base ' + \
                     'directory for dataset \"{}\".'.format(dataset)
            raise IOError(errStr)
    elif (dataset == 'CM06'):
        if (platform.system() == 'Linux'):
            basedir = '/media/jrinker/JRinker SeaGate External/data/' + \
                        'plaine-morte/CM06/'
        elif (platform.system() == 'Windows'):
            basedir = drive + '\\data\\plaine-morte\\CM06'
        if not os.path.exists(basedir):
            errStr = 'Incorrect or unavailable base ' + \
                     'directory for dataset \"{}\".'.format(dataset)
            raise IOError(errStr)

    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return basedir


#def makemetadata(dataset):
#    """ Construct metadata table
#    """
#
#    # get base directory, check it exists
#    basedir = getBasedir(dataset)
#
#    # process NREL dataset
#    if (dataset == 'NREL'):
#        metadata = makeNRELmetadata(basedir)
#
#    else:
#        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
#        raise AttributeError(errStr)
#
#    return metadata


def listmetadata(dataset,i,list_mats):
    """ Return list of parameters all heights for element i
        in list_mats
    """
    
    import os
    

    if (dataset in ['NREL','fluela','CM06']):
        IDs      = datasetSpecs(dataset)[2]             # instrument IDs
        basedir  = getBasedir(dataset)                  # base directory
        fpath    = os.path.join(basedir,list_mats[i])   # path to mat file
        n_fields = len(metadataFields(dataset))         # list wind parameters
        
        # try to load the high-frequency structure, return arrays of NaNs if failed
        try:
            struc = scio.loadmat(fpath)
        except Exception as e:
            print('Cannot load {}'.format(fpath))
            print('  ' + str(e))
            parms    = np.empty(n_fields)
            parms[:] = np.nan
            return [parms for _ in range(IDs.size)]

        # loop through sonic heights
        h_parms = []
        for ID in IDs:
            
            # calculate parameters
            row = struc2metadata(dataset,struc,ID)

            # append to list of structure parametsr
            h_parms.append(row)
            
    else:
        errStr = 'Dataset {} not coded'.format(dataset)
        raise AttributeError(errStr)
                
    return h_parms


def struc2metadata(dataset,struc_hf,ID):
    """ Calculate metadata parameters from high-frequency .mat

        Args:
            struc_hf (dictionary): high-frequency structure
            height (int): measurement height

        Returns:
            parameters (numpy array): 1D array of metadata parameters
    """
    

    if (dataset in ['NREL','fluela','CM06']):

        # get list of metadata fieldnames
        md_fields = metadataFields(dataset)

        # intialize array of metadata parameters
        parameters = np.empty(len(md_fields))

        # calculate all fields
        outdict = calculatefield(dataset,struc_hf,ID)

        # assign fields to location in output array if dictionary is non-empty
        if outdict:
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


def calculatefield(dataset,struc_hf,ID):
    """ Save atmophseric parameters in output dictionary

    """
    import calendar, time
    
    
    # meteorological constants
    g, R, kappa = 9.81, 287, 0.41
    
    # number of time steps and time increment
    N, dt = datasetSpecs(dataset)[:2]

    if (dataset == 'NREL'):

        refht = 3                   # height [m] for ref. pressure/temp/etc
        
        # get interpolation/closest heights
        T_hts = interpolationHeights(dataset,ID,'Temperature')
        loht_T, hiht_T = T_hts
        clht_T = T_hts[np.abs(np.array(T_hts)-ID).argmin()]
        loht_WS, hiht_WS = interpolationHeights(dataset,ID,'Wind_Speed_Cup')
        loht_WD, hiht_WD = interpolationHeights(dataset,ID,'Wind_Direction')
        loht_DP, hiht_DP = interpolationHeights(dataset,ID,'Dewpt_Temp')
        
        # load all necessary time series, incrementally checking flags
        clean  = 1                          # initialize clean flag
        fields = ['Sonic_u','Sonic_v','Sonic_w','Sonic_T','Temperature',\
            'Temperature','Temperature','Wind_Speed_Cup','Wind_Speed_Cup',\
            'Wind_Direction','Wind_Direction','Precipitation',\
            'Dewpt_Temp','Dewpt_Temp','Dewpt_Temp','Pressure',\
            'Temperature']
        ts_hts = [ID, ID, ID, ID, loht_T, hiht_T, clht_T, loht_WS, hiht_WS,\
                loht_WD, hiht_WD, ID, loht_DP, hiht_DP, refht, refht,\
                refht]
        time_series = np.empty((len(fields),N))
        for i in range(len(fields)):
            
            # load time series for that field, measurement height
            field = fields[i]
            ts_ht = ts_hts[i]
            outdict  = loadtimeseries(dataset,field,ts_ht,struc_hf)

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
            fname     = struc_hf['tower'][0,0][23][0,0][0][0,0][8][0]
            rec_time  = fname2time(dataset,fname)
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
            Tbar_in_K = interpolateparameter(dataset,ID,Tbar_lo_K,Tbar_hi_K,\
                                            'Temperature')
            WSbar_in  = interpolateparameter(dataset,ID,WSbar_lo,WSbar_hi,\
                                            'Wind_Speed_Cup')
            WDbar_in  = interpolateparameter(dataset,ID,WDbar_lo,WDbar_hi,\
                                            'Wind_Direction')
            Tbar_0_K  = C2K(Tbar_0)
            if (DPbar_0 > 0): A, B = 7.5, 237.3
            else:            A, B = 9.5, 265.5
            e0        = 6.11 * 10 ** ((DPbar_0*A)/(DPbar_0 + B))
            q0        = e0 / Pbar_0
            Tv0       = (Tbar_0_K)*(1 + 0.61*q0)
            dPdz      = - (g * Pbar_0) / (R * Tv0)
            DPbar_z   = interpolateparameter(dataset,ID,DPbar_lo,
                                             DPbar_hi,'Dewpt_Temp')
            Tbar_z_K  = interpolateparameter(dataset,ID,Tbar_lo_K,
                                             Tbar_hi_K,'Temperature')

            if (DPbar_z > 0): A, B = 7.5, 237.3
            else:            A, B = 9.5, 265.5
            ez        = 6.11 * 10 ** ((DPbar_z*A)/(DPbar_z + B))
            Pz        = Pbar_0 + (ID - refht)*dPdz
            qz        = ez / Pz
            Tvbar_z_K = (Tbar_z_K)*(1 + 0.61*qz)

            # initialize output dictionary
            outdict = {}

            # save values
            outdict['Record_Time']     = rec_time
            outdict['Processed_Time']  = calendar.timegm(time.gmtime())   
            outdict['Height']          = ID
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
        
        # no need to interpolate -- only have sonic measurements
        clean  = 1                          # initialize clean flag
        
        # load all necessary time series, incrementally checking flags
        fields = ['Sonic_u','Sonic_v','Sonic_w','Sonic_T',\
            'Sonic_Cup','Sonic_Direction']
        ts_hts = [ID, ID, ID, ID, ID, ID]
        time_series = np.empty((len(fields),N))
        for i in range(len(fields)):
            
            # load time series for that field, measurement height
            field = fields[i]
            ts_ht = ts_hts[i]
            outdict  = loadtimeseries(dataset,field,ts_ht,struc_hf)
            
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
            fname     = struc_hf['tower'][0,0][23][0,0][0][0,0][3][0]
            rec_time  = fname2time(dataset,fname)
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
            Tvbar_K   = Tbar_K
    
            # initialize output dictionary
            outdict = {}
    
            # save values
            outdict['Record_Time']     = rec_time
            outdict['Processed_Time']  = calendar.timegm(time.gmtime())   
            outdict['Height']          = ID
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
                                         /(kappa * g * wpTp_bar)
            outdict['MO_Length_virt']  = - (Tvbar_K * ustar**3)/ \
                                         (kappa * g * wpTp_bar)
            
        else:
            outdict = {}
            
    elif (dataset == 'CM06'):
        
        # parameters for humidity calculations
        P_1atm = 101325                     # pressure in Pascals at 1 atm
        R_dryair = 287.058                  # specfici gas constant dry air
        
        # no need to interpolate -- only have sonic measurements
        clean       = 1                     # initialize clean flag
        Pcln_thresh = 0.95                  # percentage clean for "healthy"
        
        # load all necessary time series, incrementally checking flags
        fields = ['Sonic_x','Sonic_y','Sonic_z','Sonic_T',
                  'Humidity','Humidity','Grp_Warning']
        ts_IDs = [ID, ID, ID, ID, 1, 2, ID]
        time_series = np.empty((len(fields),N))
        for i in range(len(fields)):
            
            # load time series for that field, measurement height
            field   = fields[i]
            ts_ID   = ts_IDs[i]
            outdict = loadtimeseries(dataset,field,ts_ID,struc_hf)
                        
            # save time series if it no flags
            if (len(outdict['flags']) == 0):
                time_series[i,:] = outdict['clean']
            else:
                clean = 0
                
        # check grp_warning 
        perc_clean = np.sum(time_series[ \
                            fields.index('Grp_Warning'),:] == 0.)/float(N)
        if (perc_clean < Pcln_thresh):
            clean = 0
                
        # if all time series are clean
        if clean:
                        
            # put all time series in variables
            t     = np.arange(N)*dt                 # time vector
            ux    = time_series[0,:]                # Sonic_x
            uy    = time_series[1,:]                # Sonic_y
            uz    = time_series[2,:]                # Sonic_z
            T_s   = time_series[3,:]                # Sonic_T
            hum1  = time_series[4,:]                # humidity sensor 1
            hum2  = time_series[5,:]                # humidity sensor 2
            
            # rotate sonic to u, v, w
            u, v, w = RotateTimeSeries(ux,uy,uz)
            
            # get variables necessary for later calculations
            fname     = struc_hf['name'][0]
            rec_time  = fname2time(dataset,fname)
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
            WSbar     = np.nanmean(np.sqrt(ux**2 + uy**2))
            WDbar     = np.angle(np.nanmean(np.exp( \
                                    1j*np.arctan2(uy,ux))),deg=1) % 360
            hum_bar   = np.nanmean(0.5*hum1 + 0.5*hum2)
            rho_air   = P_1atm/(R_dryair * Tbar_K)
            q         = hum_bar / rho_air / 1000
            Tvbar_K   = Tbar_K * (1 + 0.51*q)
    
            # initialize output dictionary
            outdict = {}
    
            # save values
            outdict['Record_Time']       = rec_time
            outdict['Processed_Time']    = calendar.timegm(time.gmtime())   
            outdict['ID']                = ID
            outdict['Sonic_Cup']         = WSbar
            outdict['Sonic_Direction']   = WDbar
            outdict['Mean_Wind_Speed']   = np.nanmean(u)
            outdict['Sigma_u']           = np.nanstd(up)
            outdict['Concentration_u']   = rhou
            outdict['Location_u']        = muu
            outdict['Sigma_v']           = np.nanstd(vp)
            outdict['Concentration_v']   = rhov
            outdict['Location_v']        = muv
            outdict['Sigma_w']           = np.nanstd(wp)
            outdict['Concentration_w']   = rhow
            outdict['Location_w']        = muw
            outdict['up_wp']             = upwp_bar          
            outdict['vp_wp']             = vpwp_bar
            outdict['wp_Tp']             = wpTp_bar
            outdict['up_vp']             = upvp_bar
            outdict['Tau_u']             = calculateKaimal(up + \
                                                np.nanmean(u),dt)
            outdict['Tau_v']             = calculateKaimal(vp,dt)
            outdict['Tau_w']             = calculateKaimal(wp,dt)
            outdict['Tbar_K']            = Tbar_K
            outdict['Tvbar_K']           = Tvbar_K
            outdict['MO_Length']         = -(Tbar_K * ustar**3) \
                                         /(kappa * g * wpTp_bar)
            outdict['MO_Length_virt']    = - (Tvbar_K * ustar**3)/ \
                                         (kappa * g * wpTp_bar)
            outdict['Specific_Humidity'] = q
            
        else:
            outdict = {}
        
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return outdict


def dcgain(dataset,field,ID):
    """ DC gain for specified dataset, field, and instrument ID

        Args:
            dataset (string): flag to indicate which dataset to analyze
            field (string): key for instrument
            ID (int): instrument identifier

        Returns:
            dc_gain (float): DC gain
            
    """
    
    if (dataset == 'CM06'):
        # dc_gains from "CM zeros.xls" in "Preparation and Miscellaneous Files"
        # in CM 2006 folder form Dr. Elie Bou-Zeid
        
        dc_fields = ['Sonic_x','Sonic_y','Sonic_z','Sonic_T']
        dc_gains  = np.array([[0.0210,	0.0023,	0.0045,	-0.2986],
                    [0.0071,	0.0088,	-0.0065,	0.3890],
                    [0.0298,	0.1033,	0.0210,	0.3197],
                    [0.0586,	-0.0376,	-0.0196,	0.0484],
                    [0.0689,	0.0607,	0.0173,	-0.7343],
                    [0.0134,	0.0164,	0.0033,	-0.0125],
                    [0.0646,	0.0632,	0.0184,	0.2777],
                    [0.0307,	0.0438,	-0.0109,	0.0050],
                    [0.0018,	0.0686,	0.0091,	0.4837],
                    [0.0112,	-0.0209,	-0.0125,	-0.8009],
                    [0.0458,	-0.0082,	-0.0027,	0.0579],
                    [0.1691,	0.0608,	-0.0463,	0.2648]])
        dc_gain = dc_gains[ID-1,dc_fields.index(field)]

    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return dc_gain


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
    
    
def RotateTimeSeries(ux,uy,uz):
    """ Yaw and pitch time series so v- and w-directions have zero mean
    
        Args:
            ux (numpy array): array of x-sonic velocity
            uy (numpy array): array of y-sonic velocity
            uz (numpy array): array of z-sonic velocity
            
        Returns:
            x_rot (numpy array): [n_t x 3] array of rotated data (yaw+pitch)
            x_yaw (numpy array): [n_t x 3] array of rotated data (yaw)
    """
        
    # combine velocities into array
    x_raw = np.concatenate((ux.reshape(ux.size,1),
                            uy.reshape(ux.size,1),
                            uz.reshape(ux.size,1)),axis=1)
    
    # interpolate out any NaN values
    for i_comp in range(x_raw.shape[1]):
        x               = x_raw[:,i_comp]
        idcs_all        = np.arange(x.size)
        idcs_notnan     = np.logical_not(np.isnan(x))
        x_raw[:,i_comp] = np.interp(idcs_all,
                            idcs_all[idcs_notnan],x[idcs_notnan])
        
    # rotate through yaw angle
    hyp = np.sqrt(np.mean(x_raw[:,0])**2 + np.mean(x_raw[:,1])**2)
    cos_theta = np.mean(x_raw[:,0]) / hyp
    sin_theta = np.mean(x_raw[:,1]) / hyp
    A_yaw = np.array([[cos_theta, -sin_theta, 0],
                      [sin_theta, cos_theta, 0],
                      [0, 0, 1]])
    x_yaw = np.dot(x_raw,A_yaw)
    
    # rotate through pitch angle
    hyp = np.sqrt(np.mean(x_yaw[:,0])**2 + np.mean(x_yaw[:,2])**2)
    cos_phi = np.mean(x_yaw[:,0]) / hyp
    sin_phi = np.mean(x_yaw[:,2]) / hyp
    A_pitch = np.array([[cos_phi, 0, -sin_phi],
                        [0, 1, 0],
                        [sin_phi, 0, cos_phi]])
    x_rot = np.dot(x_yaw,A_pitch)
    
    # calculate yaw angle
    yaw_theta = np.arctan2(sin_theta,cos_theta)
    
    # define rotated velocities
    u = x_rot[:,0]
    v = x_rot[:,1]
    w = x_rot[:,2]
    
#    return x_rot, x_yaw, yaw_theta
    return u,v,w
    

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


def wrappedCauchySample(shape,rho,mu):
    """ Random numbers from a Wrapped Cauchy distribution with 
        location parameter mu and concentration parameter rho.
        
        Reference: Statistical Analysis of Circular Data, Fisher, Sec. 3.3.4
        Modified to correctly implement arccos.
        
        Args:
            shape (tuple): shape of output sample
            rho (float): concentration parameter
            mu (float): location parameter
    
        Returns:
            theta (numpy array): sample of angles
    """
    
    
    if ((rho < 0) or (rho>1)):
        print 'rho must be between 0 and 1'
        return []
    
    U = np.random.rand(*shape)                      # uniform random numbers
    
    V = np.cos(2.*np.pi*U)                          # See lines before Eq. 3.28
    c = 2.*rho/(1+(rho**2))                         # See lines before Eq. 3.28
    
    B = 2*np.round(np.random.rand(*shape)) - 1     # boolean RV
                
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

        dtheta = wrappedCauchySample((n_f-1,n_m),\
            rho,mu)                                 # phase differences
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


def IEC_VelProfile(z,ZRef,URef):
    """ IEC velocity profile for values in array z with reference height
        ZRef and and reference velocitry URef
        
        Args:
            z (numpy array): array of heights in meters for velocity calcs
            ZRef (float): reference height in meters
            URef (float): reference-height mean velocity in m/s
        
        Returns:
            V (numpy array): array of mean wind speeds at heights in z
    """
    

    alpha = 0.2;

    V = URef * np.power( \
        z/ZRef, alpha );

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

      
    Iref = IEC_Iref(turbc)                      # reference hub-height TI  
    sigma1 = IEC_Sigma1(Iref,Vhub)              # reference standard deviation
    V = IEC_VelProfile(z,zhub,Vhub)             # mean wind profile
    Ti = sigma1 / V                             # TI profile
    
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
    n_f     = uniqueComponents(n_t)          # no. unique frequencies
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
    

    # read file
    tsout = io.readModel( fname )

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
    f, Cohij = vectorSpatCoh(Xi_all,Xj_all,df);

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


def WriteTurbSimInputs(fname_inp,TSDict,wr_dir):
    """ Write .inp and .spc (if UserSpec) files given TurbSim parameters
    
        Args:
            fname_inp (string): name of desired .inp file with extension
            TSDict (dictionary): dictionary with TurbSim parameters
            wr_dir (string): directory to write files to
    """
    
    import os
    
    # template filename
    tmpl_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'templates')
    spctemp = os.path.join(tmpl_dir,'Template_UsrSpc.spc')
    TStemp  = os.path.join(tmpl_dir,'Template_TurbSim.inp')
    
    # convert items in dictionary to variables for coding convenience
    T,dt,R1,R2         = TSDict['T'],TSDict['dt'],TSDict['R1'],TSDict['R2']
    DY,n_y,DZ,n_z      = TSDict['DY'],TSDict['n_y'],TSDict['DZ'],TSDict['n_z']
    URef,ZHub,ZRef     = TSDict['URef'],TSDict['ZHub'],TSDict['ZRef']
    TurbModel,ProfType = TSDict['TurbModel'],TSDict['ProfType']
    rho_u,rho_v,rho_w  = TSDict['rho_u'],TSDict['rho_v'],TSDict['rho_w']
    mu_u,mu_v,mu_w     = TSDict['mu_u'],TSDict['mu_v'],TSDict['mu_w']
    
    if (TurbModel == 'USRINP'):
        sig_u,sig_v,sig_w  = TSDict['sig_u'],TSDict['sig_v'],TSDict['sig_w']
        L_u,L_v,L_w        = TSDict['L_u'],TSDict['L_v'],TSDict['L_w']
        fname_spc          = TSDict['fname_spc']
        IECTurbc           = 'unused'
    elif (TurbModel == 'IECKAI'):
        IECTurbc           = TSDict['IECTurbc']
        fpath_spc          = 'unused'
    else:
        ValueError('No turbulence model option \"{:s}\"'.format(TurbModel))
        

    # create spectral input file if user input
    if (TurbModel == 'USRINP'):
        
        # create unscaled PSDs
        fo, ff, df = 0, 20, 0.0005                          # high-frequency spectra
        fs  = np.arange(fo,ff,df)                           # high-frequency vector
        Su    = KaimalSpectrum(fs,L_u/URef,sig_u)        # unscaled u-PSD
        Sv    = KaimalSpectrum(fs,L_v/URef,sig_v)        # unscaled v-PSD
        Sw    = KaimalSpectrum(fs,L_w/URef,sig_w)        # unscaled w-PSD
        NumF  = fs.size                                     # number of frequencies
        
        # frequenices/spectral values for scaling
        df_s, n_t = 1./T, T/dt                              # simulation parameters
        fs_s = np.arange(uniqueComponents(n_t))*df_s     # simulation freqs
        Su_s  = KaimalSpectrum(fs_s,L_u/URef,sig_u)      # simulation u-PSD
        Sv_s  = KaimalSpectrum(fs_s,L_v/URef,sig_v)      # simulation v-PSD
        Sw_s  = KaimalSpectrum(fs_s,L_w/URef,sig_w)      # simulation w-PSD
        Suk_s, Svk_s, Swk_s = Su_s*df_s,\
                        Sv_s*df_s, Sw_s*df_s                # continuous -> discrete
        alpha1 = spectralScale(Suk_s,sig_u,n_t)**2       # u scale factor
        alpha2 = spectralScale(Svk_s,sig_v,n_t)**2       # v scale factor
        alpha3 = spectralScale(Swk_s,sig_w,n_t)**2       # w scale factor
#        Scale1,Scale2,Scale3 = alpha1,alpha2,alpha3         # set scale factors
        Scale1,Scale2,Scale3 = 1,1,1         # set scale factors
        
        fpath_spc = os.path.join(wr_dir,fname_spc)
        with open(fpath_spc,'w') as f_out:
            with open(spctemp,'r') as f_temp:
                
                # print header information
                f_out.write(f_temp.readline())
                f_out.write(f_temp.readline().format(URef,sig_u,sig_v,sig_w,
                           L_u,L_v,L_w,rho_u,rho_v,rho_w,mu_u,mu_v,mu_w ))
                f_out.write(f_temp.readline())
                f_out.write(f_temp.readline().format(NumF))
                f_out.write(f_temp.readline().format(Scale1))
                f_out.write(f_temp.readline().format(Scale2))
                f_out.write(f_temp.readline().format(Scale3))
                f_out.write(f_temp.readline())
                f_out.write(f_temp.readline())
                f_out.write(f_temp.readline())
                f_out.write(f_temp.readline())
                
                # print spectra
                for i_f in range(NumF):
                    f_out.write('{:>6.4f}  '.format(fs[i_f]) + \
                        '{:>16.6f}  '.format(alpha1*Su[i_f]) + \
                        '{:>16.6f}  '.format(alpha2*Sv[i_f]) + \
                        '{:>14.6f}\n'.format(alpha3*Sw[i_f]))
    
    
    # create TurbSim input file
    fpath_inp = os.path.join(wr_dir,fname_inp)
    with open(fpath_inp,'w') as f_out:
        with open(TStemp,'r') as f_temp: 
    
            i_line = 0
            for line in f_temp:
                if i_line == 4:
                    f_out.write(line.format(R1))
                elif i_line == 5:
                    f_out.write(line.format(R2))
                elif i_line == 18:
                    f_out.write(line.format(n_z))
                elif i_line == 19:
                    f_out.write(line.format(n_y))
                elif i_line == 20:
                    f_out.write(line.format(dt))
                elif i_line == 21:
                    f_out.write(line.format(T))
                elif i_line == 23:
                    f_out.write(line.format(ZHub))
                elif i_line == 24:
                    f_out.write(line.format(DZ))
                elif i_line == 25:
                    f_out.write(line.format(DY))
                elif i_line == 30:
                    f_out.write(line.format('\"'+TurbModel+'\"'))
                elif i_line == 31:
                    f_out.write(line.format('\"'+fname_spc+'\"'))
                elif i_line == 33:
                    f_out.write(line.format('\"'+IECTurbc+'\"'))
                elif i_line == 36:
                    f_out.write(line.format('\"'+ProfType+'\"'))
                elif i_line == 38:
                    f_out.write(line.format(ZRef))
                elif i_line == 39:
                    f_out.write(line.format(URef))
                elif i_line == 63:
                    f_out.write(line.format(rho_u,mu_u))
                elif i_line == 64:
                    f_out.write(line.format(rho_v,mu_v))
                elif i_line == 65:
                    f_out.write(line.format(rho_w,mu_w))
                else:
                    f_out.write(line)
                i_line += 1
    
    return


# %%===========================================================================
# FAST ANALYSIS
# =============================================================================

def CreateTurbineDictionary(turb_name,turb_dir,BModes=1,TModes=1):
    """ Convert information in text files to dictionary containing all of the
        turbine parameters necessary to create all FAST files for a WindPACT 
        turbine
        
        Args:
            turb_name (string): turbine name
            turb_dir (string): path to turbine directory
            BModes (Boolean): whether to get blade modes from 
                              Modes v22 output, opt.
            TModes (Boolean): whether to get tower modes from 
                              Modes v22 output, opt.
            
        Returns:
            TurbDict (dictionary): turbine parameters in dictionary
    """
    import os
    

    print('\nWriting turbine dictionary' + \
        ' for {:s} to {:s}'.format(turb_name,turb_dir))
    
    # ===================== Initialize turbine dictionary =====================
    TurbDict = {}
    TurbDict['TurbName'] = turb_name
    
    # ========================= Add blade properties ==========================
    
    # load rotor data from text file
    fBldName = os.path.join(turb_dir,'parameters\\'+turb_name+'_Blades.txt')
    with open(fBldName,'r') as f:
        key,i_line = '',0
        while (key != 'DistBldProps'):
            row = f.readline().strip('\n').rstrip().split()
            key = row[0]
            if (key == 'DistBldProps'):
                values = row[1:]
            else:
                values = [float(x) for x in row[1:]]
                
            # if list is single element, pull out that element
            if (len(values) == 1):
                values = values[0]
                
            # save extracted value
            TurbDict[key] = values
            i_line += 1
            
    # calculate blade schedule from remaining table
    n_skip = i_line
    BldSched = np.genfromtxt(fBldName,skip_header=n_skip).tolist()
        
    # load mode shape data from Modes output if available
    fModesName = os.path.join(turb_dir,'modes\\'+turb_name+'_BldModes.mod')
    BldModes = np.empty((3,5))
    BldModes[:] = np.nan
    if BModes:
        with open(fModesName,'r') as f:
            i_line = 0
            for line in f:
                if ((i_line >= 20) and (i_line <= 24)):
                    BldModes[0,i_line-20] = float(line.split()[1])
                    BldModes[1,i_line-20] = float(line.split()[2])
                elif ((i_line >= 33) and (i_line <= 37)):
                    BldModes[2,i_line-33] = float(line.split()[1])
                i_line += 1
    BldModes = BldModes.tolist()
    
    # save values in dictionary
    TurbDict['BldFile']        = turb_name + '_Blade.dat'
    TurbDict['BldModes']       = BldModes
    TurbDict['BldSched']       = BldSched
        
    # ==================== Add AeroDyn properties ==================
    
    # load AeroDyn data from text file
    fADName = os.path.join(turb_dir,'parameters\\'+turb_name+'_AeroDyn.txt')
    with open(fADName,'r') as f:
        key,i_line = '',0
        while (key != 'ADSchedFields'):
            row = f.readline().strip('\n').rstrip().split()
            key = row[0]
            if (key == 'FoilNm'):
                values = [s.split('.')[0] for s in row[1:]]
            elif (key == 'ADSchedFields'):
                values = [s.split('.')[0] for s in row[1:]]
            else:
                values = [float(x) for x in row[1:]]
                
            # if list is single element, pull out that element
            if (len(values) == 1):
                values = values[0]
                
            # save extracted value
            TurbDict[key] = values
            i_line += 1
            
    # calculate AeroDyn schedule from remaining table
    n_skip             = i_line
    ADSched            = np.genfromtxt(fADName,skip_header=n_skip).tolist()
    
    # save variables
    TurbDict['ADFile']  = turb_name + '_AD.dat'
    TurbDict['ADSched'] = ADSched
        
    # ======================== Add nacelle properties =========================
    
    # load nacelle data from text file
    fNacName = os.path.join(turb_dir,'parameters\\'+turb_name+'_Nacelle.txt')
    with open(fNacName,'r') as f:
        for line in f:
            row = line.strip('\n').rstrip().split()
            if row:
                key          = row[0]
                value        = float(row[1])
                TurbDict[key] = value
            
    # ========================== Add tower properties =========================
    
    # load nacelle data from text file
    fTwrName = os.path.join(turb_dir,'parameters\\'+turb_name+'_Tower.txt')
    with open(fTwrName,'r') as f:
        key,i_line = '',0
        while (key != 'TwrSchedFields'):
            row = f.readline().strip('\n').rstrip().split()
            key = row[0]
            if (key == 'TwrSchedFields'):
                values = [s.split('.')[0] for s in row[1:]]
            else:
                values = [float(x) for x in row[1:]]
                
            # if list is single element, pull out that element
            if (len(values) == 1):
                values = values[0]
                
            # save extracted value
            TurbDict[key] = values
            i_line += 1
            
    # calculate tower schedule from remaining table
    n_skip                 = i_line
    TwrSched               = np.genfromtxt(fTwrName,skip_header=n_skip).tolist()
    
    # save variables
    TurbDict['TwrFile']  = turb_name + '_Tower.dat'
    TurbDict['TwrSched']   = TwrSched
        
    # ======================== Add control properties =========================
    
    # load control data from text file
    fCntrName = os.path.join(turb_dir,'parameters\\'+turb_name+'_Control.txt')
    with open(fCntrName,'r') as f:
        for line in f:
            row = line.strip('\n').rstrip().split()
            if row:
                key          = row[0]
                value        = float(row[1])
                TurbDict[key] = value
                
    return TurbDict

def InterpolateRotorParams(TurbDict): 
    """ Interpolate blade structural properties and aerodynamic properties
    
        Args:
            TurbDict (dictionary): turbine properties
            
        Returns:
            BldInterp (numpy array): interpolated blade properties
            ADInterp (numpy array): interpolated aerodynamic properties
    """
    
    
    # extract information from turbine dictionary
    BldSched = np.array(TurbDict['Rotor']['BldSched'])
    Fields   = TurbDict['Rotor']['BldSchedFields']
    BldEdges = np.array(TurbDict['Rotor']['BldEdges'])
    ADEdges  = np.array(TurbDict['Rotor']['ADEdges'])
    RotDiam  = TurbDict['Rotor']['RotDiam']
    HubDiam  = TurbDict['Rotor']['HubDiam']

    # place structural values in array
    GenAxLoc = BldSched[:,Fields.index('Gen. Axis Loc.')]
    StrcTwst = BldSched[:,Fields.index('Twist')]
    BMassDen = BldSched[:,Fields.index('Unit Weight')]
    FlpStff  = BldSched[:,Fields.index('EIFlap')]
    EdgStff  = BldSched[:,Fields.index('EIEdge')]
    GJStff   = BldSched[:,Fields.index('GJ')]
    EAStff   = BldSched[:,Fields.index('EA')] 
    Chord    = BldSched[:,Fields.index('Chord')] 
    AFID     = BldSched[:,Fields.index('Airfoil ID')] 
    
    # calculate intermediate parameters
    Station  = np.array(BldSched[:,Fields.index('Station')])
    BlFract  = (Station - Station[0])/(Station[-1] - Station[0])
    ACOff    = Chord*(0.25 - GenAxLoc)
    BladeLen = RotDiam/2. - HubDiam/2.
    ENodes = 0.5*(ADEdges[:-1] + ADEdges[1:])
    RNodes = ENodes*BladeLen + HubDiam/2.
    
    # interpolate blade structural parameters
    ChordInt = np.interp(BldEdges,BlFract,Chord)
    AeroCentInt = np.interp(BldEdges,BlFract,ACOff)/ChordInt + 0.25
    StrcTwstInt = np.interp(BldEdges,BlFract,StrcTwst)
    BMassDenInt = np.interp(BldEdges,BlFract,BMassDen)
    FlpStffInt  = np.interp(BldEdges,BlFract,FlpStff)
    EdgStffInt  = np.interp(BldEdges,BlFract,EdgStff)
    GJStffInt   = np.interp(BldEdges,BlFract,GJStff)
    EAStffInt   = np.interp(BldEdges,BlFract,EAStff)
    
    # place parameters in array
    BldInterp = np.empty((BldEdges.size,8))
    BldInterp[:,0] = BldEdges
    BldInterp[:,1] = AeroCentInt
    BldInterp[:,2] = StrcTwstInt
    BldInterp[:,3] = BMassDenInt
    BldInterp[:,4] = FlpStffInt
    BldInterp[:,5] = EdgStffInt
    BldInterp[:,6] = GJStffInt
    BldInterp[:,7] = EAStffInt
    
    # interpolate aerodynamic parameters
    TwstInt  = np.interp(ENodes,BlFract,StrcTwst)
    ChordInt = np.interp(ENodes,BlFract,Chord)
    AFInt    = np.round(np.interp(ENodes,BlFract,AFID))
    AFInt[1] = 2
    
    # place parameters in array
    ADInterp = np.empty((ADEdges.size-1,5))
    ADInterp[:,0] = RNodes
    ADInterp[:,1] = TwstInt
    ADInterp[:,2] = BladeLen/(ADEdges.size - 1)
    ADInterp[:,3] = ChordInt
    ADInterp[:,4] = AFInt
    
    return BldInterp, ADInterp


def InterpolateTowerParams(TurbDict): 
    """ Interpolate tower structural properties
    
        Args:
            TurbDict (dictionary): turbine properties
            
        Returns:
            TowerInterp (numpy array): interpolated tpwer properties
    """
    
    
    # extract information from turbine dictionary
    TopDiam    = TurbDict['Tower']['TopDiam']
    TopThick   = TurbDict['Tower']['TopThick'] 
    BaseDiam   = TurbDict['Tower']['BaseDiam']  
    BaseThick  = TurbDict['Tower']['BaseThick']
    ParaMass   = TurbDict['Tower']['ParaMass']
    TowerDens  = TurbDict['Tower']['TowerDens']
    TowerE     = TurbDict['Tower']['TowerE']
    TowerG     = TurbDict['Tower']['TowerG']
    TowerEdges = TurbDict['Tower']['TowerEdges']
        
    # interpolate parameters
    Diam      = np.interp(TowerEdges,[0,1],
                           [BaseDiam,TopDiam])
    Thickness = np.interp(TowerEdges,[0,1],
                           [BaseThick,TopThick])
    
    # calculate intermediate parameters
    Area         = np.pi * Diam * Thickness
    TowerMassDen = Area * TowerDens * (1.0 + ParaMass)
    EA           = Area * TowerE
    Izz          = Area * (Diam **2) / 4.
    GJ           = Izz * TowerG
    Iyy          = Izz/2
    EI           = Iyy * TowerE
    
    # place parameters in array
    TowerInterp = np.empty((len(TowerEdges),6))
    TowerInterp[:,0] =  TowerEdges
    TowerInterp[:,1] =  TowerMassDen
    TowerInterp[:,2] =  EI
    TowerInterp[:,3] =  EI
    TowerInterp[:,4] =  GJ
    TowerInterp[:,5] =  EA

    return TowerInterp


def writeBldModes(fpath_temp,fpath_out,TurbDict):
    """ Blade input file for Modes v22
    """
    
    
    # get interpolated blade structural parameters
    BldInterp = InterpolateRotorParams(TurbDict)[0]
    
    # calculate modes-specific vales
    RotorRad = TurbDict['Rotor']['RotDiam']/2.
    SSAngVel = TurbDict['Nacelle']['RatedTipSpeed']/2/np.pi/RotorRad*60.
    HubRad   = TurbDict['Rotor']['HubDiam']/2
    
    # open template file and file to write to
    with open(fpath_temp,'r') as f_temp:
        with open(fpath_out,'w') as f_write:
            i_line = 0
            for line in f_temp:
                if i_line == 1:
                    f_write.write(line.format(SSAngVel))
                elif i_line == 3:
                    f_write.write(line.format(RotorRad))
                elif i_line == 4:
                    f_write.write(line.format(HubRad))
                elif i_line == 5:
                    f_write.write(line.format(0.0))
                elif i_line == 8:
                    f_write.write(line.format(len(BldInterp)))
                elif i_line == 12:
                    for i_BlNode in range(len(BldInterp)):
                        row = [BldInterp[i_BlNode,0],0.,
                               BldInterp[i_BlNode,3],
                               BldInterp[i_BlNode,4],BldInterp[i_BlNode,5]]
                        f_write.write(line.format(*row))
                else:
                    f_write.write(line)
                i_line += 1
    
    return
    

def writeTwrModes(fpath_temp,fpath_out,TurbDict):
    """ Tower input file for Modes v22
    """
    
    # calulate interpolated tower structural properties
    TowerInterp   = InterpolateTowerParams(TurbDict)
    
    # load/calculate tower-modes-specific vales
    MainframeMass = TurbDict['Nacelle']['MainFrameMass']
    GenShaftMass  = TurbDict['Nacelle']['GenShaftMass']
    RotShaftMass  = TurbDict['Nacelle']['RotShaftMass']
    HubMass       = TurbDict['Nacelle']['HubMass']
    HubHeight     = TurbDict['Nacelle']['HubHeight']
    HHtoTT        = TurbDict['Tower']['HHtoTop']
    TowerHeight   = HubHeight - HHtoTT
    TowerTopMass  = MainframeMass + GenShaftMass + RotShaftMass + HubMass
    
    # open template file and file to write to
    with open(fpath_temp,'r') as f_temp:
        with open(fpath_out,'w') as f_write:
            i_line = 0
            for line in f_temp:
                if i_line == 3:
                    f_write.write(line.format(TowerHeight))
                elif i_line == 5:
                    f_write.write(line.format(TowerTopMass))
                elif i_line == 8:
                    f_write.write(line.format(len(TowerInterp)))
                elif i_line == 12:
                    for i_BlNode in range(len(TowerInterp)):
                        row = [TowerInterp[i_BlNode,0],
                               TowerInterp[i_BlNode,1],
                               TowerInterp[i_BlNode,2],TowerInterp[i_BlNode,3]]
                        f_write.write(line.format(*row))
                else:
                    f_write.write(line)
                i_line += 1
    
    return

def writeBlade(fpath_temp,fpath_out,TurbDict):
    """ Blade input file for FAST v7.02
    """
    
    
    # calculate blade-specific vales
    TName     = TurbDict['TName']
    title_str = 'FAST v7.02 blade file for turbine \"{:s}\" (JRinker, Duke University)'.format(TName)
    Damping   = TurbDict['Rotor']['Damping']
    RotModes  = np.array(TurbDict['Rotor']['RotModes'])
    BldInterp = InterpolateRotorParams(TurbDict)[0]
    
    # open template file and file to write to
    with open(fpath_temp,'r') as f_temp:
        with open(fpath_out,'w') as f_write:
            i_line = 0
            for line in f_temp:
                if i_line == 2:
                    f_write.write(line.format(title_str))
                elif i_line == 4:
                    f_write.write(line.format(len(BldInterp)))
                elif i_line == 6:
                    f_write.write(line.format(Damping[0]))
                elif i_line == 7:
                    f_write.write(line.format(Damping[1]))
                elif i_line == 8:
                    f_write.write(line.format(Damping[2]))
                elif i_line == 18:
                    for i_BlNode in range(len(BldInterp)):
                        f_write.write(line.format(*BldInterp[i_BlNode,:]))
                elif ((i_line >= 20) and (i_line <= 34)):
                    f_write.write(line.format(
                        RotModes.reshape(RotModes.size)[i_line-20]))
                else:
                    f_write.write(line)
                i_line += 1
    
    return
    
    
def writeTower(fpath_temp,fpath_out,TurbDict):
    """ Tower input file for FAST v7.02
    """
    
    
    # calculate tower-specific vales
    TName       = TurbDict['TName']
    title_str   = 'FAST v7.02 tower file for turbine \"{:s}\" (JRinker, Duke University)'.format(TName)
    Damping     = TurbDict['Tower']['TowerDamping']
    TwrModes    = np.array(TurbDict['Tower']['TwrModes'])
    
    # interpolate tower structural properties
    TowerInterp = InterpolateTowerParams(TurbDict)
    
    # open template file and file to write to
    with open(fpath_temp,'r') as f_temp:
        with open(fpath_out,'w') as f_write:
            i_line = 0
            for line in f_temp:
                if i_line == 2:
                    f_write.write(line.format(title_str))
                elif i_line == 4:
                    f_write.write(line.format(len(TowerInterp)))
                elif i_line == 6:
                    f_write.write(line.format(Damping[0]))
                elif i_line == 7:
                    f_write.write(line.format(Damping[1]))
                elif i_line == 8:
                    f_write.write(line.format(Damping[2]))
                elif i_line == 9:
                    f_write.write(line.format(Damping[3]))
                elif i_line == 21:
                    for i_TwrNode in range(len(TowerInterp)):
                        f_write.write(line.format(*TowerInterp[i_TwrNode,:]))
                elif ((i_line >= 23) and (i_line <= 32)):
                    f_write.write(line.format(
                        TwrModes.reshape(TwrModes.size)[i_line-23]))
                elif ((i_line >= 34) and (i_line <= 43)):
                    f_write.write(line.format(
                        TwrModes.reshape(TwrModes.size)[i_line-34]))
                else:
                    f_write.write(line)
                i_line += 1  
    
    return
    
    
def writeAeroDynTemplate(fpath_temp,fpath_out,TurbDict):
    """ AeroDyn input file for FAST v7.02
    """
    
    # calculate file-specific vales
    TurbName       = TurbDict['TurbName']
    Comment1       = 'FAST v7.02 AeroDyn file for turbine ' + \
                        '\"{:s}\" (JRinker, Duke University)'.format(TurbName)
    WindHH         = TurbDict['TowerHt'] + TurbDict['Twr2Shft'] + \
                        TurbDict['OverHang']*np.sin(TurbDict['ShftTilt']*np.pi/180.)
    
    # open template file and file to write to
    with open(fpath_temp,'r') as f_temp:
        with open(fpath_out,'w') as f_write:
            i_line = 0
            
            # read each line in template file
            for r_line in f_temp:
                
                # print comment lines
                if (i_line == 0):
                    w_line = r_line.format(Comment1)
                    f_write.write(w_line)
                    
                # print key to line if number identifier present
                elif ('{:' in r_line):
                    field = r_line.split()[1]
                    
                    # try to load value from turbine dictionary
                    try:
                        w_line = r_line.format(TurbDict[field])
                        
                    # if key isn't present, check a few conditions
                    except KeyError:
                        if (field == 'HH'):
                            w_line = r_line.format(WindHH)
                            f_write.write(w_line)
                        elif (field == 'NumFoil'):
                            # *** ENDED HERE
                            NumFoil = len(TurbDict['FoilNm'])
                            w_line = r_line.format(NumFoil)
                            f_write.write(w_line)
                            r_line = f.readlin()
#                            for i_foil in range()
                        else:
                            print(field)
                            w_line = r_line
                            f_write.write(w_line)
                    
                # copy all other lines without modification
                else:
                    w_line = r_line
                    f_write.write(w_line)
                i_line += 1  
    
    return
    
def WriteFASTTemplate(fpath_temp,fpath_out,TurbDict):
    """ FAST input file for FAST v7.02
    """
    
    # get list of dictionary fields
    DictFields = TurbDict.keys()
    
    # create header strings
    turb_name     = TurbDict['TurbName']
    Comment1 = 'FAST v7.02 input file for turbine \"{:s}\"'.format(turb_name)
    Comment2 = 'Generated by Jennifer Rinker (Duke University) ' + \
                'based on WindPACT Excel input files'.format(turb_name)
        
    # open template file and file to write to
    with open(fpath_temp,'r') as f_temp:
        with open(fpath_out,'w') as f_write:
            i_line = 0
            
            # read each line in template file
            for r_line in f_temp:
                
                # print comment lines
                if (i_line == 2):
                    w_line = r_line.format(Comment1)
                    f_write.write(w_line)
                elif (i_line == 3):
                    w_line = r_line.format(Comment2)
                    f_write.write(w_line)
                    
                # print key to line if number identifier present
                elif ('{:' in r_line):
                    field = r_line.split()[1]
                    
                    # try to load value from turbine dictionary
                    try:
                        w_line = r_line.format(TurbDict[field])
                        
                    # if key isn't present, see if a subkey is or it's a file
                    except KeyError:
                        subkey = [sk for sk in DictFields if sk in field] 
                        if subkey:
                            w_line = r_line.format(TurbDict[subkey[0]])
                        elif ('File' in field):
                            w_line = r_line.format('unused')
                        else:
                            w_line = r_line
                    f_write.write(w_line)
                    
                # copy all other lines without modification
                else:
                    w_line = r_line
                    f_write.write(w_line)
                i_line += 1   
    
    return

def writePitch(fpath_temp,fpath_out,TurbDict):
    """ Pitch controller input file for UserVSControl by ACH (FAST v7.02)
    """
    
        
    # open template file and file to write to
    with open(fpath_temp,'r') as f_temp:
        with open(fpath_out,'w') as f_write:
            i_line = 0
            
            # read each line in template file
            for r_line in f_temp:
                
                # print key to line if number identifier present
                if ('{:' in r_line):
                    field = r_line.split()[1]
                    
                    # try to load value from turbine dictionary
                    try:
                        w_line = r_line.format(TurbDict[field])
                        
                    # if key isn't present, see if a subkey is or it's a file
                    except KeyError:
                        w_line = r_line
                    f_write.write(w_line)
                    
                # copy all other lines without modification
                else:
                    w_line = r_line
                    f_write.write(w_line)
                i_line += 1   
    
    return
            
        
def writeDISCON(fpath_temp,fpath_out,TurbDict):
    """ DISCON Bladed controller routine
    """
    
    
    # load/calculate needed values
    MinPitchAng = TurbDict['Control']['MinPitchAng']
    MaxPitchAng = TurbDict['Control']['MaxPitchAng']
    KP = TurbDict['Control']['KP']
    KI = TurbDict['Control']['KI']
    KD = TurbDict['Control']['KD']
    CornerFreq = TurbDict['Control']['CornerFreq']
    MaxPitchRate = TurbDict['Control']['MaxPitchRate']
    GSStartAng = TurbDict['Control']['GSStartAng']
    GSEndAng = TurbDict['Control']['GSEndAng']
    GSExp = TurbDict['Control']['GSExp']
    PCTimeStep = TurbDict['Control']['PCTimeStep']
    VSTimeStep = TurbDict['Control']['VSTimeStep']
    CutInGenSpd = TurbDict['Control']['CutInGenSpd']
    MaxTrqRate = TurbDict['Control']['MaxTrqRate']
    GenTrqAlpha = TurbDict['Control']['GenTrqAlpha']
    MaxTrq = TurbDict['Control']['MaxTrq']
    Rgn15Spd = TurbDict['Control']['Rgn15Spd']
    Rgn3MinPitch = TurbDict['Control']['Rgn3MinPitch']
    Rgn2SlpPrc = TurbDict['Control']['Rgn2SlpPrc']
    RatedGenRPM = TurbDict['Nacelle']['RatedGenRPM']

    # get intermediate values
    PC_RefSpd = RatedGenRPM*np.pi/30
    GSCoef    = 1./(GSStartAng**GSExp)
    RatedPower = TurbDict['Nacelle']['RatedPower']
    GenEff     = TurbDict['Nacelle']['GenEff']
                
   # open template file and file to write to
    with open(fpath_temp,'r') as f_temp:
        with open(fpath_out,'w') as f_write:
            i_line = 0
            for line in f_temp:
                if i_line == 31:
                    f_write.write(line.format(CornerFreq))
                elif i_line == 44:
                    f_write.write(line.format(PCTimeStep))
                elif i_line == 45:
                    f_write.write(line.format(GSStartAng))
                elif i_line == 46:
                    f_write.write(line.format(GSEndAng))
                elif i_line == 47:
                    f_write.write(line.format(GSCoef))
                elif i_line == 48:
                    f_write.write(line.format(GSExp))
                elif i_line == 49:
                    f_write.write(line.format(KD))
                elif i_line == 50:
                    f_write.write(line.format(KI))
                elif i_line == 51:
                    f_write.write(line.format(KP))
                elif i_line == 52:
                    f_write.write(line.format(MaxPitchAng))
                elif i_line == 53:
                    f_write.write(line.format(MaxPitchRate))
                elif i_line == 54:
                    f_write.write(line.format(MinPitchAng))
                elif i_line == 55:
                    f_write.write(line.format(PC_RefSpd))
                elif i_line == 66:
                    f_write.write(line.format(CutInGenSpd))
                elif i_line == 67:
                    f_write.write(line.format(VSTimeStep))
                elif i_line == 68:
                    f_write.write(line.format(MaxTrqRate))
                elif i_line == 69:
                    f_write.write(line.format(MaxTrq))
                elif i_line == 70:
                    f_write.write(line.format(GenTrqAlpha))
                elif i_line == 71:
                    f_write.write(line.format(Rgn15Spd))
                elif i_line == 72:
                    f_write.write(line.format(Rgn3MinPitch*np.pi/180))
                elif i_line == 73:
                    f_write.write(line.format(RatedGenRPM*np.pi/30))
                elif i_line == 74:
                    f_write.write(line.format(RatedPower/GenEff))
                elif i_line == 77:
                    f_write.write(line.format(Rgn2SlpPrc))
                else:
                    f_write.write(line)
                i_line += 1
    
    return
    
    
def writeFASTFiles(turb_dir,TName,wind_fname,
                   BlPitch0=None,RotSpeed0=None,
                   wind_dir=None,fileID=None,TMax=630.0,
                   wr_dir=None):
    """ Copy FAST template files from subfolder ``templates'' in turb_dir and
        place into turb_dir, pasting in nessecary initial conditions and
        simulation values as necessary.
        
        Args:
            turb_dir (string): path to turbine FAST directory
            TName (string): turbine name
            wind_fname (string): name of wind file
            wind_dir (string): optional path to directory with wind files
            BlPitch0 (list/numpy array): optional initial blade pitch angles
            RotSpee0 (float): optional initial rotor speed
            fileID (string): optional file identifier
            TMax (float): maximum simulation time
            
    """
    import os
    
    
    
    # get initial wind speed
    wind_fpath = os.path.join(wind_dir,wind_fname)
    u0 = GetFirstWind(wind_fpath)
    
    print('Writing FAST files for \"{:s}\" '.format(TName) + \
            'with wind file {:s}'.format(wind_fpath))
    
    # set optional values as necessary
    GenDOF = 'True'
    if wind_dir is None:
        wind_dir = os.path.join(turb_dir,'Wind')
    if wr_dir is None:
        wr_dir = turb_dir
    if fileID is None:
        fAD_name  = wind_fname[:-4] + '_AD.ipt'
        fFST_name = wind_fname[:-4] + '.fst'
    else:
        fAD_name  = TName+'_'+fileID+'_AD.ipt'
        fFST_name = TName+'_'+fileID+'.fst'
    if BlPitch0 is None:
        mdict = scio.loadmat(os.path.join(turb_dir,'steady_state',
                                        TName+'_SS.mat'),squeeze_me=True)
        LUT        = mdict['SS']
        saveFields = [str(s).strip() for s in mdict['Fields']]
        BlPitch0 = np.interp(u0,LUT[:,saveFields.index('WindVxi')],
                            LUT[:,saveFields.index('BldPitch1')])*np.ones(3)
    if RotSpeed0 is None:
        mdict = scio.loadmat(os.path.join(turb_dir,'steady_state',
                                        TName+'_SS.mat'),squeeze_me=True)
        LUT        = mdict['SS']
        saveFields = [str(s).strip() for s in mdict['Fields']]
        RotSpeed0 = np.interp(u0,LUT[:,saveFields.index('WindVxi')],
                            LUT[:,saveFields.index('RotSpeed')])

    # create filenames
    fAD_temp  = os.path.join(turb_dir,'templates',TName+'_AD.ipt')
    fAD_out   = os.path.join(wr_dir,fAD_name)
    fFST_temp = os.path.join(turb_dir,'templates',TName+'.fst')
    fFST_out  = os.path.join(wr_dir,fFST_name)
    
    # write AeroDyn file
    with open(fAD_temp,'r') as f_temp:
        with open(fAD_out,'w') as f_write:
            i_line = 0
            for line in f_temp:
                if i_line == 9:
                    f_write.write(line.format(wind_fpath))
                else:
                    f_write.write(line)
                i_line += 1
                
    # write FAST file
    with open(fFST_temp,'r') as f_temp:
        with open(fFST_out,'w') as f_write:
            i_line = 0
            for line in f_temp:
                if i_line == 9:
                    f_write.write(line.format(TMax))
                elif i_line == 45:
                    f_write.write(line.format(BlPitch0[0]))
                elif i_line == 46:
                    f_write.write(line.format(BlPitch0[1]))
                elif i_line == 47:
                    f_write.write(line.format(BlPitch0[2]))
                elif i_line == 59:
                    f_write.write(line.format(GenDOF))
                elif i_line == 72:
                    f_write.write(line.format(RotSpeed0))
                elif i_line == 160:
                    f_write.write(line.format(fAD_name))
                else:
                    f_write.write(line)
                i_line += 1
                
    return
    
    
def PlotSSTurbineResponse(x,Data,Fields,fig=None):
    """ Steady-state turbine response in 3-axis plot
    
        Args:
            x (numpy array): mean wind speed
            Data (numpy array): array of turbine data
            Fields (numpy array): columns in Data
            SS (boolean): whether response is steady-state or time series
            fig (pyplot figure handle): handle to figure
            
        Returns:
            fig (pyplot figure handle): handle to figure
    """
    import matplotlib.pyplot as plt
    
    
    # create figure if none specified
    if fig is None:
        fig = plt.figure(figsize=(7,9))
        
    # axes locations
    xPlot = 0.11
    yPlot = np.arange(3)[::-1]*0.31 + 0.08
    wd,ht = 0.80,0.25
    
    # Axes 1: GenSpeed,RotPwr,GenPwr,RotThrust,RotTorq
    ax1 = fig.add_axes([xPlot,yPlot[0],wd,ht])
    
    plt.plot(x,Data[:,Fields.index('GenSpeed')],label='GenSpeed, rpm')
    plt.plot(x,Data[:,Fields.index('LSShftPwr')],label='RotPwr, kW')
    plt.plot(x,Data[:,Fields.index('GenPwr')],label='GenPwr, kW')
    plt.plot(x,Data[:,Fields.index('RotThrust')],label='RotThrust, kN')
    plt.plot(x,Data[:,Fields.index('RotTorq')],label='RotTorq, kN-m')
    
    ax1.set_xlim([x[0],x[-1]])
    ax1.grid('on')
    ax1.legend(loc=2,fontsize='small')
    ax1.set_xticks(np.arange(x[0],x[-1]+1))
    
    
    # Axes 2: RotSpeed,BlPitch,GenTq,TSR
    ax2 = fig.add_axes([xPlot,yPlot[1],wd,ht])
    
    plt.plot(x,Data[:,Fields.index('RotSpeed')],label='RotSpeed, rpm')
    plt.plot(x,Data[:,Fields.index('BldPitch1')],label='BlPitch, $^\mathrm{o}$')
    plt.plot(x,Data[:,Fields.index('GenTq')],label='GenTq, kN-m')
    plt.plot(x,Data[:,Fields.index('TSR')],label='TSR, -')
    
    ax2.set_xlim([x[0],x[-1]])
    ax2.grid('on')
    ax2.legend(loc=2,fontsize='small')
    ax2.set_xticks(np.arange(x[0],x[-1]+1))
    
    # Axes 3: OopDefl1,IPDefl1,TTDspFA,TTDspSS
    ax3 = fig.add_axes([xPlot,yPlot[2],wd,ht])
    
    plt.plot(x,Data[:,Fields.index('OoPDefl1')],label='OoPDefl1, m')
    plt.plot(x,Data[:,Fields.index('IPDefl1')],label='IPDefl1, m')
    plt.plot(x,Data[:,Fields.index('TTDspFA')],label='TTDspFA, m')
    plt.plot(x,Data[:,Fields.index('TTDspSS')],label='TTDspSS, m')
    
    #ax3.set_ylim([0,40])
    ax3.set_xlim([x[0],x[-1]])
    ax3.grid('on')
    ax3.legend(loc=2,fontsize='small')
    ax3.set_xticks(np.arange(x[0],x[-1]+1))
    
    return fig,ax1,ax2,ax3
    
    
def PlotTurbineResponse(t,Data,Fields,fig=None):
    """ Turbine response in 3-axis plot
    
        Args:
            t (numpy array): time
            Data (numpy array): array of turbine data
            Fields (numpy array): columns in Data
            SS (boolean): whether response is steady-state or time series
            fig (pyplot figure handle): handle to figure
            
        Returns:
            fig (pyplot figure handle): handle to figure
    """
    import matplotlib.pyplot as plt
    
    
    # create figure if none specified
    if fig is None:
        fig = plt.figure(figsize=(6.5,9))
    fig.clf()
    
    # axes locations
    xPlot = 0.10
    yPlot = [0.88,0.586,0.32,0.06]
    wd,ht = 0.83,0.17
    
    # Axes 0: WindVxi
    ax0 = fig.add_axes([xPlot,yPlot[0],wd,0.06])
    
    plt.plot(t,Data[:,Fields.index('WindVxi')])
    
    ax0.set_xlim([t[0],t[-1]])
    ax0.set_xticklabels([])
    ax0.set_ylabel('Wind [m/s]')
    plt.locator_params(axis='y',nbins=4)
    
    # Axes 1: GenSpeed,RotPwr,GenPwr,RotThrust,RotTorq
    ax1 = fig.add_axes([xPlot,yPlot[1],wd,ht])
    
    plt.plot(t,Data[:,Fields.index('GenSpeed')],label='GenSpeed, rpm')
    plt.plot(t,Data[:,Fields.index('LSShftPwr')],label='RotPwr, kW')
    plt.plot(t,Data[:,Fields.index('GenPwr')],label='GenPwr, kW')
    plt.plot(t,Data[:,Fields.index('RotThrust')],label='RotThrust, kN')
    plt.plot(t,Data[:,Fields.index('RotTorq')],label='RotTorq, kN-m')
    
    ax1.set_xlim([t[0],t[-1]])
    ax1.legend(fontsize='small',
               bbox_to_anchor=(1.0, 1.02), loc=4, borderaxespad=0.)
#    ax1.set_xlabel('Time [s]')
    ax1.set_xticklabels([])
    
    
    # Axes 2: RotSpeed,BlPitch,GenTq,TSR
    ax2 = fig.add_axes([xPlot,yPlot[2],wd,ht])
    
    plt.plot(t,Data[:,Fields.index('RotSpeed')],label='RotSpeed, rpm')
    plt.plot(t,Data[:,Fields.index('BldPitch1')],label='BlPitch, $^\mathrm{o}$')
    plt.plot(t,Data[:,Fields.index('GenTq')],label='GenTq, kN-m')
#    plt.plot(t,Data[:,Fields.index('TSR')],label='TSR, -')
    
    ax2.set_xlim([t[0],t[-1]])
    ax2.legend(fontsize='small',
               bbox_to_anchor=(1.0, 1.02), loc=4, borderaxespad=0.)
#    ax2.set_xlabel('Time [s]')
    ax2.set_xticklabels([])
    
    # Axes 3: OopDefl1,IPDefl1,TTDspFA,TTDspSS
    ax3 = fig.add_axes([xPlot,yPlot[3],wd,ht])
    
    plt.plot(t,Data[:,Fields.index('OoPDefl1')],label='OoPDefl1, m')
    plt.plot(t,Data[:,Fields.index('OoPDefl2')],label='OoPDefl2, m')
    plt.plot(t,Data[:,Fields.index('OoPDefl3')],label='OoPDefl3, m')
#    plt.plot(t,Data[:,Fields.index('OoPDefl1')],label='OoPDefl32, m')
#    plt.plot(t,Data[:,Fields.index('IPDefl1')],label='IPDefl1, m')
#    plt.plot(t,Data[:,Fields.index('TTDspFA')],label='TTDspFA, m')
#    plt.plot(t,Data[:,Fields.index('TTDspSS')],label='TTDspSS, m')
    
    #ax3.set_ylim([0,40])
    ax3.set_xlim([t[0],t[-1]])
    ax3.legend(fontsize='small',
               bbox_to_anchor=(1.0, 1.02), loc=4, borderaxespad=0.)
    ax3.set_xlabel('Time [s]')
    
    return fig,ax1,ax2,ax3
    
    
def writeSteadyWind(U,wind_dir=''):
    """ Write steady-state wind file
    
        Args:
            U (float): wind value
            wind_dir (string): optional path to directory with wind files
    """
    import os
    
    # create path to file
    wind_fname = 'NoShr_'+'{:2.1f}'.format(U).zfill(4)+'.wnd'
    wind_fpath = os.path.join(wind_dir,wind_fname)
    
    # write wind file
    with open(wind_fpath,'w') as f:
        f.write('! Wind file for steady {:.1f} m/s wind without shear.\n'.format(U))
        f.write('! Time	Wind	Wind	Vert.	Horiz.	Vert.	LinV	Gust\n')
        f.write('!	Speed	Dir	Speed	Shear	Shear	Shear	Speed\n')
        f.write('   0.0\t{:.1f}\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n'.format(U))
        f.write('   0.1\t{:.1f}\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n'.format(U))
        f.write('9999.9\t{:.1f}\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n'.format(U))
        
    return
    
    
def writeStepWind(U0,UF,t_step=40.0,dt_step=0.1,wind_dir=''):
    """ Write step wind file
    
        Args:
            U0 (float): starting wind value
            UF (float): ending wind value
            t_step (float): optional time step occurs
            dt_step (float): optional time step duration
            wind_dir (string): optional path to directory with wind files
    """
    import os
    
    # create path to file
    U0s, UFs   = '{:2.1f}'.format(U0).zfill(4), '{:2.1f}'.format(UF).zfill(4)
    wind_fname = 'Step_'+U0s+'_'+UFs+'.wnd'
    wind_fpath = os.path.join(wind_dir,wind_fname)
    
    # write wind file
    with open(wind_fpath,'w') as f:
        f.write('! Wind file for step wind ({:.1f}'.format(U0) + \
                ' - {:.1f} m/s) without shear.\n'.format(UF))
        f.write('! Time	Wind	Wind	Vert.	Horiz.	Vert.	LinV	Gust\n')
        f.write('!	Speed	Dir	Speed	Shear	Shear	Shear	Speed\n')
        f.write('   0.0\t{:.1f}\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n'.format(U0))
        f.write('   0.1\t{:.1f}\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n'.format(U0))
        f.write('{:6.1f}\t{:.1f}\t0.0\t0.0\t0.0'.format(t_step,U0) + \
                    '\t0.0\t0.0\t0.0\n')
        f.write('{:6.1f}\t{:.1f}\t0.0\t0.0\t0.0'.format(t_step+dt_step,UF) + \
                    '\t0.0\t0.0\t0.0\n')
        f.write('9999.9\t{:.1f}\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n'.format(UF))
        
    return
    
    
def writeKaimalWind(U,sig,tau,rho,mu,T=600.0,dt=0.05,T_steady=30.0,
                    fileID='',wind_dir=''):
    """ 1D Kaimal wind file (single point)
    
        Args:
            U (float): mean wind speed
            sigma (float): turbulent standard deviation
            tau (float): Kaimal time scale
            rho (float): concentration parameter
            mu (float): location parameter
            T (float): optional length of turbulent portion of simulation
            dt (float): optional time step
            T_steady (float): optional length of steady time before turbulent
            fileID (string): optional unique file identifier
            wind_dir (string): optional path to directory with wind files
    """
    import os
    
    
    # simulate wind, add time before
    n_t = np.ceil(T/dt)
    t,u = generateKaimal1D(n_t,1,dt,U,sig,tau,rho,mu)
    t_all = np.arange((T+T_steady)/dt)*dt
    u_all = np.empty(t_all.shape)
    u_all[:T_steady/dt] = u[0]
    u_all[T_steady/dt:] = np.squeeze(u)    
    
    # create path to file
    if len(fileID)>0:
        wind_fname = 'Kaimal'+'_'+fileID+'.wnd'
    else:
        wind_fname = 'Kaimal.wnd'
    wind_fpath = os.path.join(wind_dir,wind_fname)
    
    # write wind file
    with open(wind_fpath,'w') as f:
        f.write('! Wind file for 1D Kaimal simulation: ' + \
                '{:.2f} {:.2f} {:.2f} {:.2f} {:.2f}\n'.format(U,sig,tau,rho,mu))
        f.write('! Time	Wind	Wind	Vert.	Horiz.	Vert.	LinV	Gust\n')
        f.write('!	Speed	Dir	Speed	Shear	Shear	Shear	Speed\n')
        for i_t in range(t_all.size):
            f.write('{:6.2f}\t{:.2f}\t0.0\t'.format(t_all[i_t],u_all[i_t]) + \
                    '0.0\t0.0\t0.0\t0.0\t0.0\n')

    return
    
    
def writeHarmonicWind(U,A,freq,T=630.0,dt=0.05,T_steady=30.0,
                    n_osc=None,fileID='',wind_dir=''):
    """ 1D Kaimal wind file (single point)
    
        Args:
            U (float): mean wind speed
            sigma (float): turbulent standard deviation
            tau (float): Kaimal time scale
            rho (float): concentration parameter
            mu (float): location parameter
            T (float): optional length of turbulent portion of simulation
            dt (float): optional time step
            T_steady (float): optional length of steady time before turbulent
            fileID (string): optional unique file identifier
            wind_dir (string): optional path to directory with wind files
    """
    import os
    
    
    # initialize time, wind speed vectors
    t = np.arange(0,T+dt,dt)
    u = U*np.ones(t.shape)
    i_steady = int(T_steady/dt)
    if n_osc is None:
        u[i_steady:] = A*np.sin(2*np.pi*freq*(t[i_steady:]-T_steady))
    else:
        i_end = i_steady + int(n_osc*(1./freq)/dt)
        u[i_steady:i_end] = U + A*np.sin(2*np.pi*freq* \
                                            (t[i_steady:i_end]-T_steady))
        
    # create path to file
    if len(fileID)>0:
        wind_fname = 'Harmonic'+'_'+fileID+'.wnd'
    else:
        wind_fname = 'Harmonic.wnd'
    wind_fpath = os.path.join(wind_dir,wind_fname)
    
    # write wind file
    with open(wind_fpath,'w') as f:
        f.write('! Wind file for 1D harmonic simulation: ' + \
                '{:.2f} {:.2f} {:.2f}\n'.format(U,A,freq))
        f.write('! Time	Wind	Wind	Vert.	Horiz.	Vert.	LinV	Gust\n')
        f.write('!	Speed	Dir	Speed	Shear	Shear	Shear	Speed\n')
        for i_t in range(t.size):
            f.write('{:6.2f}\t{:.2f}\t0.0\t'.format(t[i_t],u[i_t]) + \
                    '0.0\t0.0\t0.0\t0.0\t0.0\n')
    return
    
    
def GetFirstWind(wind_fpath):
    """ First wind speed from file
    
        Args:
            wind_fpath (string): path to wind file
            
        Returns:
            u0 (float): initial wind value
    """
    import pyts.io.main as io

    # if it's a .wnd file (text)
    if wind_fpath.endswith('.wnd'):
        with open(wind_fpath,'r') as f:
            for i in range(3):
                f.readline()
            u0 = float(f.readline().split()[1])
    
    # if it's a .bts file
    elif wind_fpath.endswith('.bts'):
        
        tsout = io.readModel(wind_fpath)            # read file
        u0    = tsout.u[:,:,0].mean()               # average first wind

    else:
        errStr = 'Can only analyze .wnd (text) and .bts files'
        ValueError(errStr)

    return u0
    
    
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


def fname2time(dataset,fname):
    """ Convert filename to float representing time stamp

        Args:
            dataset (string): flag to indicate dataset
            fname (string): name of file

        Returns:
            timestamp (float): float representing timestamp
    """
    
    if (dataset == 'NREL'):

        # extract date information
        year     = int(fname[6:10])
        month    = int(fname[0:2])
        day      = int(fname[3:5])
        hour     = int(fname[11:13])
        minute   = int(fname[14:16])
        time_tup = (year,month,day,hour,minute)
    
        # convert tuple to float
        time_flt = timetup2flt(time_tup)
        
    elif (dataset == 'CM06'):

        # extract date information
        year     = int(fname[6:10])
        month    = int(fname[0:2])
        day      = int(fname[3:5])
        hour     = int(fname[11:13])
        minute   = int(fname[13:15])
        time_tup = (year,month,day,hour,minute)
    
        # convert tuple to float
        time_flt = timetup2flt(time_tup)
        
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return time_flt


def time2fpath(dataset,timestamp):
    """ Convert float to path to data structure

        Args:
            timestamp (float or tuple): float or tuple (year,month,day,hour,minute)
                                        representing record start time

        Returns:
            fpath (string): path to data 
            
    """
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


def field2datfield(dataset,field,ID):
    """ Dataset-specific fieldname corresponding with custom fieldname
        and height

        Args:
            dataset (string): flag for dataset
            field (string): custom fieldname
            ID (float): instrument identifier (often msmnt height)

        Returns:
            datfield (string): dataset-specific fieldname
    """

    ID = int(ID)                        # ID must be int for string calculations

    if dataset == 'NREL':

        if   (field == 'Sonic_u'):        datfield = 'Sonic_u_' + str(ID) + 'm'
        elif (field == 'Sonic_v'):        datfield = 'Sonic_v_' + str(ID) + 'm'
        elif (field == 'Sonic_w'):        datfield = 'Sonic_w_' + str(ID) + 'm'
        elif (field == 'Sonic_T'):        datfield = 'Sonic_Temp_rotated_' + str(ID) + 'm'
        elif (field == 'Temperature'):    datfield = 'Air_Temp_' + str(ID) + 'm'
        elif (field == 'Wind_Speed_Cup'): datfield = 'Cup_WS_' + str(ID) + 'm'
        elif (field == 'Wind_Direction'): datfield = 'Vane_WD_' + str(ID) + 'm'
        elif (field == 'Precipitation'):  datfield = 'PRECIP_INTEN'
        elif (field == 'Dewpt_Temp'):     datfield = 'Dewpt_Temp_' + str(ID) + 'm'
        elif (field == 'Pressure'):       datfield = 'Baro_Presr_3m'
        else:
            errStr = 'Unknown custom field {} for dataset \"{}\"'.format(field,dataset)
            raise AttributeError(errStr)
            
        # check that datafield is valid
        if (not check_datfields(dataset,datfield)):
            raise ValueError('Invalid height {} for field {}'.format(\
                ID,field))

    elif dataset == 'fluela':

        if   (field == 'Sonic_u'):         datfield = 'Sonic_u_' + str(ID) + 'm'
        elif (field == 'Sonic_v'):         datfield = 'Sonic_v_' + str(ID) + 'm'
        elif (field == 'Sonic_w'):         datfield = 'Sonic_w_' + str(ID) + 'm'
        elif (field == 'Sonic_T'):         datfield = 'Sonic_Temp_rotated_' + str(ID) + 'm'
        elif (field == 'Sonic_Cup'):       datfield = 'Sonic_CupEqHorizSpeed_' \
                                                         + str(ID) + 'm'
        elif (field == 'Sonic_Direction'): datfield = 'Sonic_direction_' \
                                                         + str(ID) + 'm'
        else:
            errStr = 'Unknown custom field {} for dataset \"{}\"'.format(field,dataset)
            raise AttributeError(errStr)
            
        # check that datafield is valid
        if (not check_datfields(dataset,datfield)):
            raise ValueError('Invalid height {} for field {}'\
                             .format(ID,field))
                             
    elif dataset == 'CM06':
        
        grp_nmbr,inst_nmbr = id2grpnmbr(dataset,ID)

        if   (field == 'Sonic_x'):     datfield = 'Ux_{:d}_{:d}(1,1)'.format(grp_nmbr,inst_nmbr)
        elif (field == 'Sonic_y'):     datfield = 'Uy_{:d}_{:d}(1,1)'.format(grp_nmbr,inst_nmbr)
        elif (field == 'Sonic_z'):     datfield = 'Uz_{:d}_{:d}(1,1)'.format(grp_nmbr,inst_nmbr)
        elif (field == 'Sonic_T'):     datfield = 'Ts_{:d}_{:d}(1,1)'.format(grp_nmbr,inst_nmbr)
        elif (field == 'Humidity'):    datfield = 'kh2o_{:d}_1'.format(ID)
        elif (field == 'Grp_Warning'): datfield = 'grp_warning_{:d}(1)'.format(grp_nmbr)
        else:
            errStr = 'Unknown custom field {} for dataset \"{}\"'.format(field,dataset)
            raise AttributeError(errStr)
            
        # check that datafield is valid
        if (not check_datfields(dataset,datfield)):
            raise ValueError('Invalid identifier {:d} for field {:s}'\
                             .format(ID,field))

    else:
        errStr = 'Dataset \"{}\" not coded.'.format(dataset)
        raise KeyError(errStr)

    return datfield


def id2grpnmbr(dataset,ID):
    """ Instrument ID to group and instrument number
    
        Args:
            dataset (string): flag for dataset
            ID (int): instrument ID

        Returns:
            grp_nmbr (int): group number
            inst_nmbr (int): instrument number
    """
    
    if (dataset == 'CM06'):
        grp_nmbr  = (ID-1) / 4 + 1
        inst_nmbr = (ID-1) % 4 + 1
        
    else:
        errStr = 'Dataset \"{}\" not coded.'.format(dataset)
        raise KeyError(errStr)
    
    return grp_nmbr, inst_nmbr


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
