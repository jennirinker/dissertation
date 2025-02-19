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
import os, re, glob, time, datetime, platform, json
import numpy as np
import scipy.io as scio
from peakdetect import peakdetect
from rainflow import rainflow
import statsmodels.api as sm
import scipy.stats

# =============================================================================
# ---------------------------- FILE I/O ---------------------------------------
# =============================================================================

def loadmetadata(dataset):
    """ Load data and fields from metadata file for dataset
    
        Args:
            dataset (string): flag for dataset
            
        Returns:
            fields (list): fields in each column
            metadata (numpy array): values for each field and each record
    """
    
    # hard-code base location of data
    BaseDir = 'C:\\Users\\jrinker\\Dropbox\\research\\processed_data'
    
    # define path to mat file
    MatName = '{:s}-metadata.mat'.format(dataset)
    MatPath = os.path.join(BaseDir,MatName)
    
    # load mat file
    struc     = scio.loadmat(MatPath)       # load structure
    fields    = [str(string.rstrip(' ')) for \
              string in  struc['fields']]   # load fields, stripping spaces
    metadata  = struc['metadata']
    
    # clean metadata
    metadata = cleanmetadata(dataset,metadata,fields)

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
        
    # calculate datfield and make sure it's valid before proceeding
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
        specs   = datasetSpecs(dataset)
        N, dt   = specs['n_t'], specs['dt']
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
        
    elif (dataset == 'PM06'):

        # set desired length of time series
        specs   = datasetSpecs(dataset)
        N, dt   = specs['n_t'], specs['dt']
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
        if ('Grp' not in datfield):     
            
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
    
            # remove spikes/detrend sonics, remove DC gain
            # detrend everything else
            if (len(flags) == 0):
                if any([key in datfield for key in ('Sonic_x','Sonic_y',
                                                    'Sonic_z','Sonic_T')]):
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

    elif (dataset == 'texastech'):
        
        # set data ranges, desired length of time series
        dataRng = dataRanges(dataset,datfield)
        specs   = datasetSpecs(dataset)
        N, dt   = specs['n_t'], specs['dt']
        t       = np.arange(N)*dt
                    
        # determine what kind of data was provided
        if (isinstance(data_in,dict)):                   # data structure given
            struc = data_in
        elif (isinstance(data_in,(int,float,tuple))):    # timestamp given
            fpath = time2fpath(dataset,data_in)
	    # try to load the structure
	    try:
		struc = scio.loadmat(fpath,squeeze_me=True)
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
    Units = [s.replace('\xb7','-') for s in Units]  # replace weird character
    
    # save everything in the dictionary
    FASTDict['Data']   = Data
    FASTDict['Fields'] = Fields
    FASTDict['Units']  = Units
    
    return FASTDict


def LoadFASTStats(RunName,TurbName,Stat,Parm,
                  scale=1):
    """ Load FAST statistics
    
        Args:
            RunName (string): run name
            stat (string): statistic to load
            parm (string): parameter of interes
            scale (float): inverse scaling factor for output data [opt]
            
        Returns:
            x (numpy array): [n_files x n_p] array of input data
            y (numpy array): [n_files] array of stats data
    """
    
    BaseStatDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                'processed_data\\proc_stats'
    m_strs = ['l','m','h']              # lo, middle, high for DEL
    Parms = RunName2WindParms(RunName)
    WindParms = [Parms['URefs'],Parms['Is'],
                 np.log10(Parms['Ls']),Parms['rhos']]
    
    # load the output data
    if 'DEL' in Stat:
        DictPath   = os.path.join(BaseStatDir,RunName,TurbName + '_DELstats.mat')
        stats_dict = scio.loadmat(DictPath,squeeze_me=True)
        fnames     = [s.rstrip() for s in stats_dict['fnames']]
        fields     = [s.rstrip() for s in stats_dict['DELKeys']]
        CalcStats  = stats_dict['DELStats']
        
        m_idx      = m_strs.index(Stat[-1])
        ms         = np.array([float(k.split('=')[1]) for k in fields if Parm in k])
        ParmKeys   = [k for k in fields if Parm in k]           # unsorted
        ParmKeys   = [ParmKeys[i] for i in ms.argsort()]        # sorted
        StatField  = ParmKeys[m_idx]
        
        y = CalcStats[:,fields.index(StatField)]

    else:
        DictPath   = os.path.join(BaseStatDir,RunName,TurbName + '_stats.mat')
        stats_dict = scio.loadmat(DictPath,squeeze_me=True)
        fnames     = [s.rstrip() for s in stats_dict['fnames']]
        fields     = [s.rstrip() for s in stats_dict['fields']]         # RootMOoP, etc.
        calc_stats = [s.rstrip() for s in stats_dict['calc_stats']]     # min, max, etc.
        proc_stats = stats_dict['proc_stats']                           # array of data
        
        n_fields = len(fields)
        y = proc_stats[:,calc_stats.index(Stat)*n_fields + \
                         fields.index(Parm)]                # output data
              
    # create the input data
    x = np.empty((y.size,4))
    for i_f in range(y.size):
        file_id =  fnames[i_f].rstrip('.out').split('_')[1]
        for i_p in range(len(WindParms)):
            x[i_f,i_p] = WindParms[i_p][int(file_id[i_p],16)]   # hex to int
    
    # scale output data
    y = y / float(scale)
    
    return x, y


# =============================================================================
# ------------------------ METADATA PROCESSING --------------------------------
# =============================================================================


def metadataFields(dataset):
    """ Define list of fields to be stored in metadata table

        Args:
            dataset (string): flag for dataset

        Returns:
            fields (list): list of strings defining metadata
                columns
    """
    
    # NOTE: these define the order of the metadata array, so keep order logical

    if (dataset == 'NREL'):
        fields = ['Record_Time','Processed_Time','ID','Wind_Speed_Cup', \
              'Wind_Direction','Precipitation', \
              'Mean_Wind_Speed', \
              'Sigma_u','Concentration_u','Location_u', \
              'Sigma_v','Concentration_v','Location_v', \
              'Sigma_w','Concentration_w','Location_w', \
              'up_wp','vp_wp','wp_Tp','up_vp','Tau_u','Tau_v','Tau_w', \
              'MO_Length_interp','MO_Length_near','MO_Length_virt']
    elif (dataset == 'fluela'):
        fields = ['Record_Time','Processed_Time','ID','Sonic_Cup', \
              'Sonic_Direction', 'Mean_Wind_Speed', \
              'Sigma_u','Concentration_u','Location_u', \
              'Sigma_v','Concentration_v','Location_v', \
              'Sigma_w','Concentration_w','Location_w', \
              'up_wp','vp_wp','wp_Tp','up_vp','Tau_u','Tau_v','Tau_w', \
              'MO_Length','MO_Length_virt']
    elif (dataset == 'PM06'):
        fields = ['Record_Time','Processed_Time','ID', \
              'Sonic_Cup', 'Sonic_Direction', 'Specific_Humidity',\
              'Mean_Wind_Speed', \
              'Sigma_u','Concentration_u','Location_u', \
              'Sigma_v','Concentration_v','Location_v', \
              'Sigma_w','Concentration_w','Location_w', \
              'up_wp','vp_wp','wp_Tp','up_vp','Tau_u','Tau_v','Tau_w', \
              'MO_Length','MO_Length_virt','Tbar_K','Tvbar_K']
    elif (dataset == 'texastech'):
        fields = ['Record_Time','Processed_Time','ID', \
              'Wind_Speed_UVW','Wind_Direction_UVW', \
              'Wind_Direction_Sonic','Wind_Speed_Sonic', \
              'Sigma_u_UVW','Concentration_u_UVW','Location_u_UVW', \
              'Sigma_v_UVW','Concentration_v_UVW','Location_v_UVW', \
              'Sigma_w_UVW','Concentration_w_UVW','Location_w_UVW', \
              'Tau_u_UVW','Tau_v_UVW','Tau_w_UVW', \
              'Mean_Wind_Speed', \
              'Sigma_u','Concentration_u','Location_u', \
              'Sigma_v','Concentration_v','Location_v', \
              'Sigma_w','Concentration_w','Location_w', \
              'up_wp','vp_wp','wp_Tp','up_vp','Tau_u','Tau_v','Tau_w', \
              'Specific_Humidity', 'Pressure', \
              'MO_Length','MO_Length_virt','Tbar_K','Tvbar_K']
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
    
    # NOTE: these do not need to be in any particular order

    if (dataset == 'NREL'):
        datfields = ['Sonic_u_15m','Sonic_u_30m','Sonic_u_50m','Sonic_u_76m',\
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
        datfields = ['Sonic_u_36m','Sonic_u_54m','Sonic_u_75m',\
        'Sonic_v_36m','Sonic_v_54m','Sonic_v_75m',\
        'Sonic_w_36m','Sonic_w_54m','Sonic_w_75m',\
        'Sonic_Temp_rotated_36m','Sonic_Temp_rotated_54m',\
        'Sonic_Temp_rotated_75m','Sonic_CupEqHorizSpeed_36m',\
        'Sonic_CupEqHorizSpeed_54m','Sonic_CupEqHorizSpeed_75m',\
        'Sonic_direction_36m','Sonic_direction_54m',\
        'Sonic_direction_75m']
    elif (dataset == 'PM06'):
        datfields = ['Sonic_z_12', 'Sonic_z_11', 'Sonic_z_10', 'Sonic_T_7', \
        'Sonic_T_1', 'Sonic_T_3', 'Sonic_T_2', 'Sonic_T_5', 'Sonic_T_4', \
        'Hygro_1', 'Sonic_T_6', 'Sonic_T_9', 'Sonic_T_8', 'Sonic_y_1', \
        'Sonic_z_7', 'Sonic_z_4', 'Sonic_T_11', 'Sonic_T_10', 'Sonic_T_12', \
        'Sonic_x_11', 'Sonic_x_10', 'Sonic_x_12', 'Sonic_y_4', 'Sonic_y_5', \
        'Sonic_y_6', 'Sonic_y_7', 'Sonic_x_9', 'Sonic_x_8', 'Sonic_y_2', \
        'Sonic_y_3', 'Sonic_x_5', 'Sonic_x_4', 'Sonic_x_7', 'Sonic_x_6', \
        'Sonic_x_1', 'Sonic_y_9', 'Sonic_x_3', 'Sonic_x_2', 'Sonic_z_9', \
        'Sonic_z_8', 'Time', 'Sonic_z_3', 'Sonic_z_2', 'Sonic_z_1', \
        'Sonic_y_12', 'Sonic_z_6', 'Sonic_y_10', 'Sonic_y_11', 'Sonic_y_8', \
        'Sonic_z_5', 'GrpWarn_1', 'Hygro_3', 'GrpWarn_3', 'GrpWarn_2', 'Hygro_2']
    elif (dataset == 'texastech'):
        datfields = ['Air_Temp_10m', 'Air_Temp_116m', 'Air_Temp_158m', 'Air_Temp_17m', \
        'Air_Temp_1m', 'Air_Temp_200m', 'Air_Temp_2m', 'Air_Temp_47m', \
        'Air_Temp_4m', 'Air_Temp_75m', 'Baro_Presr_10m', 'Baro_Presr_116m', \
        'Baro_Presr_158m', 'Baro_Presr_17m', 'Baro_Presr_1m', 'Baro_Presr_200m', \
        'Baro_Presr_2m', 'Baro_Presr_47m', 'Baro_Presr_4m', 'Baro_Presr_75m', \
        'Rel_Hum_10m', 'Rel_Hum_116m', 'Rel_Hum_158m', 'Rel_Hum_17m', \
        'Rel_Hum_1m', 'Rel_Hum_200m', 'Rel_Hum_2m', 'Rel_Hum_47m', \
        'Rel_Hum_4m', 'Rel_Hum_75m', 'Sonic_T_10m', 'Sonic_T_116m', \
        'Sonic_T_158m', 'Sonic_T_17m', 'Sonic_T_1m', 'Sonic_T_200m', \
        'Sonic_T_2m', 'Sonic_T_47m', 'Sonic_T_4m', 'Sonic_T_75m', \
        'Sonic_x_10m', 'Sonic_x_116m', 'Sonic_x_158m', 'Sonic_x_17m', \
        'Sonic_x_1m', 'Sonic_x_200m', 'Sonic_x_2m', 'Sonic_x_47m', 'Sonic_x_4m', \
        'Sonic_x_75m', 'Sonic_y_10m', 'Sonic_y_116m', 'Sonic_y_158m', \
        'Sonic_y_17m', 'Sonic_y_1m', 'Sonic_y_200m', 'Sonic_y_2m', \
        'Sonic_y_47m', 'Sonic_y_4m', 'Sonic_y_75m', 'Sonic_z_10m', \
        'Sonic_z_116m', 'Sonic_z_158m', 'Sonic_z_17m', 'Sonic_z_1m', \
        'Sonic_z_200m', 'Sonic_z_2m', 'Sonic_z_47m', 'Sonic_z_4m', \
        'Sonic_z_75m', 'UVW_x_10m', 'UVW_x_116m', 'UVW_x_158m', 'UVW_x_17m', \
        'UVW_x_200m', 'UVW_x_47m', 'UVW_x_4m', 'UVW_x_75m', 'UVW_y_10m', \
        'UVW_y_116m', 'UVW_y_158m', 'UVW_y_17m', 'UVW_y_200m', 'UVW_y_47m', \
        'UVW_y_4m', 'UVW_y_75m', 'UVW_z_10m', 'UVW_z_116m', 'UVW_z_158m', \
        'UVW_z_17m', 'UVW_z_200m', 'UVW_z_47m', 'UVW_z_4m', 'UVW_z_75m']

    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)
        
    return datfield in datfields


def datasetSpecs(dataset):
    """ Measurement specifics for dataset

        Args:
            dataset (string): flag for dataset

        Returns:
            n_t (int): number of time steps per record
            dt (float): time step
            IDs (numpy array): anemometer IDs (height or number)
    """
    
    specs = {}
    
    if (dataset == 'NREL'):
        specs['n_t']          = 12000
        specs['dt']           = 0.05
        specs['IDs']          = np.array([15,30,50,76,100,131])
        specs['WSlim']        = 3.0
        specs['dir1']         = 240.
        specs['dir2']         = 330.
        specs['Prelim']       = 2.7
    elif (dataset == 'fluela'):
        specs['n_t']          = 6000
        specs['dt']           = 0.10
        specs['IDs']          = np.array([36,54,75])
        specs['WSlim']        = 3.0
        specs['dir1']         = 225.
        specs['dir2']         = 165.
        specs['T1']           = timetup2flt((2010,1,1,0,0))
        specs['T2']           = timetup2flt((2010,3,23,0,0))
    elif (dataset == 'PM06'):
        specs['n_t']          = 12000
        specs['dt']           = 0.05
        specs['IDs']          = np.arange(1,13)
        specs['sonic_offset'] = 120.
        specs['WSlim']        = 3.0
        specs['dir1']         = 60.
        specs['dir2']         = 180.
    elif (dataset == 'texastech'):
        specs['n_t']          = 30000
        specs['dt']           = 0.02
        specs['IDs']          = np.array([1,2,4,10,17,47,75,116,158,200])
        specs['sonic_offset'] = 30.
        specs['UVW_offset']   = 75.
        specs['WSlim']        = 3.0
        specs['dir1']         = 210.
        specs['dir2']         = 30.
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return specs


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
        if ('Sonic_u' in datfield):     # m/s
            dataRng = [-35.05,35.05]
        elif ('Sonic_v' in datfield):   # m/s
            dataRng = [-35.05,35.05]
        elif ('Sonic_w' in datfield):   # m/s
            dataRng = [-29.95,29.95]
        elif ('Sonic_T' in datfield):   # C
            dataRng = [-19.95,49.95]
        elif ('Air_Temp' in datfield):  # C
            dataRng = [-50.,50.]
        elif ('Cup_WS' in datfield):    # m/s
            dataRng = [0.,80.]
        elif ('Vane_WD' in datfield):   # degrees
            dataRng = [0.,360.]
        elif ('PRECIP' in datfield):    # none
            dataRng = [0.,4.]
        elif ('Dewpt' in datfield):     # C
            dataRng = [-50.,50.]
        elif ('Presr' in datfield):     # mbar
            dataRng = [740.,1000.]
        else:
              raise KeyError('Field {} not recognized.'.format(datfield))
              
    elif (dataset == 'fluela'):

        # define data ranges
        if ('Sonic_u' in datfield):     # m.s
            dataRng = [-35.05,35.05]
        elif ('Sonic_v' in datfield):   # m/s
            dataRng = [-35.05,35.05]
        elif ('Sonic_w' in datfield):   # m/s
            dataRng = [-29.95,29.95]
        elif ('Sonic_T' in datfield):   # C
            dataRng = [-19.95,49.95]
        elif ('Sonic_Cup' in datfield): # m/s
            dataRng = [  0.00,49.57]            # sqrt(2*(35.05**2))
        elif ('Sonic_direction' in datfield):   # degrees
            dataRng = [ 0.00,360.00]            # angle from 0 to 360
        else:
              raise KeyError('Field {} not recognized.'.format(datfield))
              
    elif (dataset == 'PM06'):

        # define data ranges
        #    Sonic data ranges taken from CSAT 3 specifications
        #    Hygrometer data from KH20 specifications
        if ('Sonic_x' in datfield):     # m/s
            dataRng = [-29.95,29.95]
        elif ('Sonic_y' in datfield):   # m/s
            dataRng = [-29.95,29.95]
        elif ('Sonic_z' in datfield):   # m/s
            dataRng = [-7.95,7.95]
        elif ('Sonic_T' in datfield):   # C
            dataRng = [-29.95,49.95]
        elif ('Hygro' in datfield):     # g/m^3
            dataRng = [  1.75,19.25] 
        else:
              raise KeyError('Field {} not recognized.'.format(datfield))
              
    elif (dataset == 'texastech'):

        # define data ranges
        #    Sonic data ranges taken from CSAT 3 specifications
        #    Hygrometer data from KH20 specifications
        if ('Sonic_x' in datfield):     # m/s
            dataRng = [-29.95,29.95]
        elif ('Sonic_y' in datfield):   # m/s
            dataRng = [-29.95,29.95]
        elif ('Sonic_z' in datfield):   # m/s
            dataRng = [-7.95,7.95]
        elif ('UVW_x' in datfield):     # m/s
            dataRng = [-40,40]
        elif ('UVW_y' in datfield):     # m/s
            dataRng = [-40,40]
        elif ('UVW_z' in datfield):     # m/s
            dataRng = [-40,40]
        elif ('Sonic_T' in datfield):   # C
            dataRng = [-29.95,49.95]
        elif ('Rel_Hum' in datfield):
            dataRng = [0,100]               # %
        elif ('Air_Temp' in datfield):
            dataRng = [-50.,50.]            # C
        elif ('Baro_Presr' in datfield):
            dataRng = [74000.,100000.]      # Pa
        else:
              raise KeyError('Field {} not recognized.'.format(datfield))

    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return dataRng


def getBasedir(dataset):
    """ Get path to directory with high-frequency mat files and check if it 
        exists

        Args:
            dataset (string): flag for dataset
            drive (string): drive identifier [opt]

        Returns:
            basedir (str): path to top level of data directory
    """
    
    # define platform-specific drives to search      
    if (platform.system() == 'Linux'):
        drives = ['/media/jrinkerJRinker SeaGate External/',
                  '/media/jrinkerSeagate Backup Plus Drive/']
        
    elif (platform.system() == 'Windows'):
        drives = ['E:\\','G:\\','H:\\','T:\\','V:\\','D:\\']
        
    # get string list of drives for possible error message
    drives_str = ''
    for i in range(len(drives)-1): drives_str += drives[i] + ', '
    drives_str += 'or {:s}'.format(drives[-1])
                    
    # define dataset locations on drive
    if (dataset == 'NREL'):
        dataset_loc = 'data\\nrel-20Hz'
            
    elif (dataset == 'fluela'):
        dataset_loc = 'data\\fluela-high_freq\\'
            
    elif (dataset == 'PM06'):
        dataset_loc = 'data\\plaine-morte\\CM06\\'
            
    elif (dataset == 'texastech'):
        dataset_loc = 'data\\texas-tech\\'
        
    # loop through drives, searcing for dataset
    for drive in drives:
        basedir = drive + dataset_loc
        if os.path.exists(basedir):
            return basedir

    # if dataset not found on any drives, throw error
    errStr = 'Unable to locate data for dataset \"{:s}\"'.format(dataset) \
                + ' in drives {:s}'.format(drives_str)
    raise IOError(errStr)

    return


def listmetadata(dataset,i,list_mats):
    """ Return list of parameters all heights for element i
        in list_mats
    """
    
    # Andy-processed files can't have squeeze_me=True
    if (dataset in ['NREL','fluela']):
        IDs      = datasetSpecs(dataset)['IDs']         # instrument IDs
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
            
    # ones I processed need squeeze_me=True
    elif (dataset in ['PM06','texastech']):
        IDs      = datasetSpecs(dataset)['IDs']         # instrument IDs
        basedir  = getBasedir(dataset)                  # base directory
        fpath    = os.path.join(basedir,list_mats[i])   # path to mat file
        n_fields = len(metadataFields(dataset))         # list wind parameters
        
        # try to load the high-frequency structure, return arrays of NaNs if failed
        try:
            struc = scio.loadmat(fpath,squeeze_me=True)
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
    
    if (dataset in ['NREL','fluela','PM06','texastech']):

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
        
    # if time series is all NaNs, return nan
    if all(np.isnan(x)):
        return np.nan

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
    import calendar
    
    # meteorological constants
    g, R, kappa = 9.81, 287, 0.41
    
    # number of time steps and time increment
    specs = datasetSpecs(dataset)
    N, dt = specs['n_t'], specs['dt']

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
            P_0   = time_series[15,:]               # ref ht pressure [mbar]
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
            e0        = 6.11 * 10 ** ((DPbar_0*A)/(DPbar_0 + B)) # [mbar]
            q0        = 0.622 * e0 / Pbar_0
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
            qz        = 0.622 * ez / Pz
            Tvbar_z_K = (Tbar_z_K)*(1 + 0.61*qz)

            # initialize output dictionary
            outdict = {}

            # save values
            outdict['Record_Time']     = rec_time
            outdict['Processed_Time']  = calendar.timegm(time.gmtime())   
            outdict['ID']              = ID
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
            outdict['Tau_u']           = calculateKaimal(up + \
                                                np.nanmean(u),dt)
            outdict['Tau_v']           = calculateKaimal(vp,dt)
            outdict['Tau_w']           = calculateKaimal(wp,dt)
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
            outdict['ID']              = ID
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
            outdict['Tau_u']           = calculateKaimal(up + \
                                                np.nanmean(u),dt)
            outdict['Tau_v']           = calculateKaimal(vp,dt)
            outdict['Tau_w']           = calculateKaimal(wp,dt)
            outdict['MO_Length']       = -(Tbar_K * ustar**3) \
                                         /(kappa * g * wpTp_bar)
            outdict['MO_Length_virt']  = - (Tvbar_K * ustar**3)/ \
                                         (kappa * g * wpTp_bar)
            
        else:
            outdict = {}
            
    elif (dataset == 'PM06'):
        
        # parameters for humidity calculations
        P_1atm   = 101325.                  # pressure in Pascals at 1 atm
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
                print(field,ts_ID,outdict['flags'])
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
            WD_offset = datasetSpecs(dataset)['sonic_offset']
            WD        = WD_offset*np.pi/180. - np.arctan2(uy,ux)
            WDbar     = np.angle(np.nanmean(np.exp( \
                                    1j*WD)),deg=1) % 360
            hum_bar   = np.nanmean(0.5*hum1 + 0.5*hum2)
            rho_air   = P_1atm/(R_dryair * Tbar_K)
            q         = hum_bar / rho_air / 1000
            Tvbar_K   = Tbar_K * (1 + 0.61*q)
    
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

    elif (dataset == 'texastech'):
        
        # no need to interpolate - temperature and pressure msmnt at all heights
        clean  = 1                          # initialize clean flag
        
        # set fields of time series to load
        if ID > 8:
            fields = ['Sonic_x','Sonic_y','Sonic_z','Sonic_T',
                      'Humidity','Temperature','Pressure',
                      'UVW_x','UVW_y','UVW_z']
        else:
            fields = ['Sonic_x','Sonic_y','Sonic_z','Sonic_T',
                      'Humidity','Temperature','Pressure']
        ts_IDs = np.ones(len(fields)) * ID
        
        # load cleaned time series, incrementally checking flags
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
                
        # if all time series are clean
        if clean:
            
            # put all time series in variables
            t     = np.arange(N)*dt                 # time vector
            ux_s  = time_series[0,:]                # Sonic_x
            uy_s  = time_series[1,:]                # Sonic_y
            uz_s  = time_series[2,:]                # Sonic_z
            T_s   = time_series[3,:]                # Sonic_T
            RH    = time_series[4,:]                # humidity sensor
            T     = time_series[5,:]                # air temperature in C
            P     = time_series[6,:]                # barometric pressure [Pa]
            P     = P * 0.01                        # barometric pressure [mbar]
            if ID > 8:
                ux_p = time_series[7,:]            # propellor anemometer - x
                uy_p = time_series[8,:]            # propellor anemometer - y
                uz_p = time_series[9,:]            # propellor anemometer - z
            else:
                ux_p = np.zeros(ux_s.shape)
                uy_p = np.zeros(uy_s.shape)
                uz_p = np.zeros(uz_s.shape)
                ux_p[:] = np.nan
                uy_p[:] = np.nan
                uz_p[:] = np.nan
            
            # rotate sonic to u, v, w
            u_s, v_s, w_s = RotateTimeSeries(ux_s,uy_s,uz_s)
            
            # rotate UVW to u, v, w
            u_p, v_p, w_p = RotateTimeSeries(ux_p,uy_p,uz_p)
                        
            # get dictionary values
            fname     = struc_hf['name']
            rec_time  = fname2time(dataset,fname)
            T_K       = C2K(T)
            T_s_K     = C2K(T_s)
            up_p      = nandetrend(t,u_p)
            vp_p      = nandetrend(t,v_p)
            wp_p      = nandetrend(t,w_p)
            up_s      = nandetrend(t,u_s)
            vp_s      = nandetrend(t,v_s)
            wp_s      = nandetrend(t,w_s)
            Tp_s      = nandetrend(t,T_s_K)
            upwp_bar  = np.nanmean(up_s*wp_s)
            upvp_bar  = np.nanmean(up_s*vp_s)
            vpwp_bar  = np.nanmean(vp_s*wp_s)
            wpTp_bar  = np.nanmean(wp_s*Tp_s)
            ustar     = ( (upwp_bar)**2 + (vpwp_bar)**2 ) ** 0.25
            rhou_s,muu_s  = signalPhaseCoherence(up_s)
            rhov_s,muv_s  = signalPhaseCoherence(vp_s)
            rhow_s,muw_s  = signalPhaseCoherence(wp_s)
            rhou_p,muu_p  = signalPhaseCoherence(up_p)
            rhov_p,muv_p  = signalPhaseCoherence(vp_p)
            rhow_p,muw_p  = signalPhaseCoherence(wp_p)
            Tbar      = np.nanmean(T)                           # [C]
            Tbar_K    = np.nanmean(T_K)                         # [K]
            WSbar_s   = np.nanmean(np.sqrt(ux_s**2 + uy_s**2))
            WSbar_p   = np.nanmean(np.sqrt(ux_p**2 + uy_p**2))
            WD_s      = np.arctan2(uy_s,ux_s)
            WD_p      = np.arctan2(uy_p,ux_p)
            WDbar_s   = np.angle(np.nanmean(np.exp( \
                                    1j*WD_s)),deg=1)
            WDbar_p   = np.angle(np.nanmean(np.exp( \
                                    1j*WD_p)),deg=1)
            WDbar_s   = (180. + 30 - WDbar_s) % 360.        # convert to cardinal
            WDbar_p   = (180. + 75 - WDbar_p) % 360.        # convert to cardinal
            RHbar     = np.nanmean(RH)
            Pbar      = np.nanmean(P)                   # [mbar]
            if Tbar >= 0: A, B = 7.5, 237.3
            else:         A, B = 9.5, 265.5
            esbar     = 6.11 * 10**((A*Tbar) / (Tbar + B))
            ebar      = (RHbar / 100.) * esbar          # [mbar]
            q         = 0.622 * ebar / Pbar
            Tvbar_K   = Tbar_K * (1 + 0.61*q)

            # initialize output dictionary
            outdict = {}
    
            # save values
            outdict['Record_Time']          = rec_time
            outdict['Processed_Time']       = calendar.timegm(time.gmtime())   
            outdict['ID']                   = ID
            outdict['Wind_Speed_Sonic']     = WSbar_s
            outdict['Wind_Speed_UVW']       = WSbar_p
            outdict['Wind_Direction_Sonic'] = WDbar_s
            outdict['Wind_Direction_UVW']   = WDbar_p
            outdict['Sigma_u_UVW']          = np.nanstd(up_p)
            outdict['Concentration_u_UVW']  = rhou_p
            outdict['Location_u_UVW']       = muu_p
            outdict['Sigma_v_UVW']          = np.nanstd(vp_p)
            outdict['Concentration_v_UVW']  = rhov_p
            outdict['Location_v_UVW']       = muv_p
            outdict['Sigma_w_UVW']          = np.nanstd(wp_p)
            outdict['Concentration_w_UVW']  = rhow_p
            outdict['Location_w_UVW']       = muw_p
            outdict['Tau_u_UVW']            = calculateKaimal(up_p + \
                                                np.nanmean(u_p),dt)
            outdict['Tau_v_UVW']            = calculateKaimal(vp_p,dt)
            outdict['Tau_w_UVW']            = calculateKaimal(wp_p,dt)
            outdict['Mean_Wind_Speed']      = np.nanmean(u_s)
            outdict['Sigma_u']              = np.nanstd(up_s)
            outdict['Concentration_u']      = rhou_s
            outdict['Location_u']           = muu_s
            outdict['Sigma_v']              = np.nanstd(vp_s)
            outdict['Concentration_v']      = rhov_s
            outdict['Location_v']           = muv_s
            outdict['Sigma_w']              = np.nanstd(wp_s)
            outdict['Concentration_w']      = rhow_s
            outdict['Location_w']           = muw_s
            outdict['up_wp']                = upwp_bar          
            outdict['vp_wp']                = vpwp_bar
            outdict['wp_Tp']                = wpTp_bar
            outdict['up_vp']                = upvp_bar
            outdict['Tau_u']                = calculateKaimal(up_s + \
                                                np.nanmean(u_s),dt)
            outdict['Tau_v']                = calculateKaimal(vp_s,dt)
            outdict['Tau_w']                = calculateKaimal(wp_s,dt)
            outdict['Tbar_K']               = Tbar_K
            outdict['Tvbar_K']              = Tvbar_K
            outdict['Pressure']             = Pbar
            outdict['MO_Length']            = -(Tbar_K * ustar**3) \
                                         /(kappa * g * wpTp_bar)
            outdict['MO_Length_virt']       = - (Tvbar_K * ustar**3)/ \
                                         (kappa * g * wpTp_bar)
            outdict['Specific_Humidity']    = q
            
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
    
    if (dataset == 'PM06'):
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
    
    # return all NaNs if any component is all nan values
    if all(np.isnan(ux))*all(np.isnan(uy))*all(np.isnan(uz)):
        u = np.zeros(ux.shape)
        v = np.zeros(uy.shape)
        w = np.zeros(uz.shape)
        u[:] = np.nan
        v[:] = np.nan
        w[:] = np.nan
        
    # if at least one data point in all three components
    else:
        
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
    

# =============================================================================
# ----------------------- METADATA ANALYSIS -----------------------------------
# =============================================================================


def metadataFpath(dataset):
    """ Filepath to metadata table
    
        Args:
            dataset (string): toggle for dataset choice
            
        Returns:
            fpath (string): filepath to metadata file
    """
    
    
    base  =  'C:\\Users\\jrinker\\Dropbox\\research\\processed_data'
    fpath = os.path.join(base,dataset + '-metadata.mat')
        
    return fpath


def screenmetadata(fields,metadata,dataset):
    """ Screen the metadata for data quality
    
        Args:
            fields (list): list of fields associated with metadata cols
            metadata (numpy array): array of metadata
            dataset (string): toggle for dataset choice
            
        Returns:
            screendata (numpy array): numpy array with screened data
    """
    
    # get dataset-specific specifications (include screening parameters)
    specs = datasetSpecs(dataset)
    
    if (dataset == 'NREL'):
        CSLim  = specs['WSlim']             # lower cup speed limit
        dir1   = specs['dir1']              # CCW edge for direction range
        dir2   = specs['dir2']              # CW edge for direction range
        preLim = specs['Prelim']            # lower precipitation limit
        
        # column indices for each value
        CScol  = fields.index('Wind_Speed_Cup')
        dirCol = fields.index('Wind_Direction')
        preCol = fields.index('Precipitation')
        
        # screen data
        screendata = metadata[np.where(metadata[:,CScol] > CSLim)]
        screendata = screendata[np.where( (screendata[:,dirCol] - dir1) % 360. \
                                        < (dir2 - dir1) % 360.)]
        screendata = screendata[np.where(screendata[:,preCol] >= preLim)]

    elif (dataset == 'fluela'):
        CSLim  = specs['WSlim']             # lower cup speed limit
        dir1   = specs['dir1']              # CCW edge for direction range
        dir2   = specs['dir2']              # CW edge for direction range
        T1 = specs['T1']                    # start time
        T2 = specs['T2']                    # end time
        
        # column indices for each value
        CScol   = fields.index('Sonic_Cup')
        dirCol  = fields.index('Sonic_Direction')
        timeCol = fields.index('Record_Time')
        
        # screen data
        screendata = metadata[np.where(metadata[:,CScol] > CSLim)]
        screendata = screendata[np.where( (screendata[:,dirCol] - dir1) % 360. \
                                        < (dir2 - dir1) % 360.)]
        screendata = screendata[np.where(screendata[:,timeCol] >= T1)]
        screendata = screendata[np.where(screendata[:,timeCol] <= T2)]

    elif (dataset == 'PM06'):
        CSLim  = specs['WSlim']             # lower cup speed limit
        dir1   = specs['dir1']              # CCW edge for direction range
        dir2   = specs['dir2']              # CW edge for direction range
        
        # column indices for each value
        CScol   = fields.index('Sonic_Cup')
        dirCol  = fields.index('Sonic_Direction')
        
        # screen data
        screendata = metadata[np.where(metadata[:,CScol] > CSLim)]
        screendata = screendata[np.where( (screendata[:,dirCol] - dir1) % 360. \
                                        < (dir2 - dir1) % 360.)]
                                        
    elif (dataset == 'texastech'):
        CSLim  = specs['WSlim']             # lower cup speed limit
        dir1   = specs['dir1']              # CCW edge for direction range
        dir2   = specs['dir2']              # CW edge for direction range
        
        # column indices for screening values
        UsCol   = fields.index('Wind_Speed_Sonic')
        dirsCol = fields.index('Wind_Direction_Sonic')
        
        # screen data (by sonic to screen at all heights)
        screendata = metadata[np.where(metadata[:,UsCol] > CSLim)]
        screendata = screendata[(screendata[:,dirsCol] - dir1) % 360. \
                       < (dir2 - dir1) % 360.,:]
        
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)
        
    return screendata
    

def cleanmetadata(dataset,raw_parms,fields):
    """ Manually remove bad data (NaNs or other from Texas Tech dataset)
    
        Args:
            dataset (string): flag to indicate dataset
            raw_parms (numpy array): raw metadata
            fields (list): string indicators of fields
            
        Returns:
            cleandata (numpy array): cleaned metadata
    """
    
    # lots of manual cleaning if texas tech
    if dataset == 'texastech':
        
        # get dataset-specific specifications (include screening parameters)
        specs = datasetSpecs('texastech')
    
        # column indices for screening values
        UsCol   = fields.index('Wind_Speed_Sonic')
        UpCol   = fields.index('Wind_Speed_UVW')
        sigCol  = fields.index('Sigma_u')
        rhoCol  = fields.index('Concentration_u')
        dirpCol = fields.index('Wind_Direction_UVW')
        dirsCol = fields.index('Wind_Direction_Sonic')
        
        # filter out _only_ sonic data with NaNs (UVW has NaNs 1/2/4 m)
        metadata = raw_parms[np.logical_not(np.isnan(raw_parms[:,UsCol])),:]
        
        # manually screen for values that exceed thresholds
        UsMax  = np.nanmax(1.10 * metadata[:,UpCol])
        sigMax, rhoMax = 4., 0.6
        cleandata = metadata[metadata[:,UsCol] < UsMax,:]      # too-high wind speeds
        cleandata = cleandata[cleandata[:,sigCol] < sigMax,:]  # removes some rain
        cleandata = cleandata[cleandata[:,rhoCol] < rhoMax,:]
        
        # manually change wind directions
        phis, phip = specs['sonic_offset'],specs['UVW_offset']
        cleandata[:,dirsCol] = (180. + phis - cleandata[:,dirsCol]) % 360
        cleandata[:,dirpCol] = (180. + phip - cleandata[:,dirpCol]) % 360
        
    # just remove NaNs if not texas tech
    elif (dataset in ['NREL','fluela','PM06']):
        
        # filter out the rows with NaN values
        cleandata = raw_parms[np.logical_not( \
            np.any(np.isnan(raw_parms),axis=1)),:]
        
    else:
        ErrStr = 'Dataset \"{:s}\" not coded.'.format(dataset)
        raise ValueError(ErrStr)

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
    
    # turn p into 1d numpy array in case it's a float
    float_in = 0
    if not np.array(x).shape:
        x = np.array([x])
        float_in = 1
        
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
    
    # return to float if float was passed in
    if float_in: F_comp = F_comp[0]
        
    return F_comp


def compositeISF(p,dist_name,p_main,x_T=float("inf"),p_GP=(0.1,0,1)):
    """ Inverse survival function of single/composite
        distribution

        Args:
            p (numpy array): values at which to evaluate CDF
            dist_name (string): distribution type
            p_main (tuple): main distribution parameters
            x_T (float): optional, threshold value
            p_GP (tuple): optional, GP distribution parameters

        Returns:
            x (numpy array): inverse survival function at quantiles Q
    """
    
    # turn p into 1d numpy array in case it's a float
    float_in = 0
    if not np.array(p).shape:
        p = np.array([p])
        float_in = 1
    
    # initialize array for output
    x_comp = np.empty(p.shape)
    
    # initialize distributions
    dist_main = getattr(scipy.stats, dist_name)
    dist_GP   = getattr(scipy.stats, 'genpareto')
    
    # define CDF functions to improve code readability
    G_main = lambda Q: dist_main.isf(Q, *p_main[:-2], \
            loc=p_main[-2], scale=p_main[-1])
    G_GP   = lambda Q: dist_GP.isf(Q, *p_GP[:-2], \
            loc=p_GP[-2], scale=p_GP[-1])

    # calculate threshold quantile
    Q_T = dist_main.cdf(x_T, *p_main[:-2], \
            loc=p_main[-2], scale=p_main[-1])
    
    # get indices for main and GP distributions
    i_main, i_GP = np.where(p <= Q_T), np.where(p > Q_T)
    
    # set distribution values
    x_comp[i_main] = G_main(p[i_main])
    x_comp[i_GP]   = G_GP((p[i_GP]-Q_T)/(1-Q_T))
    
    # return to float if float was passed in
    if float_in: x_comp = x_comp[0]
    
    return x_comp


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
    elif (dist_name == 'gumbel_r'):
        bnds = ((None,None),(0,None))
    else:
        ErrStr = 'Distribution \"{:s}\" not coded'.format(dist_name)
        raise ValueError(ErrStr)
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
        Coh = IEC_SpatialCoherence(zhub,Vhub,DR,f)
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
    

def SampleWindParameters(NumSamps,dataset,BaseDir,ParmSample,iH,
                         URange=[-float('inf'),float('inf')]):
    """ Draw a correlated sample of wind parameters for a given dataset.
        Throw out wind speed values that are outside the given range.
    """
    
    SigRange = [0,float('inf')]
    TauRange   = [1e-1,1e4]
    RhoRange = [0,1]
    
    # load marginal distribution information
    dist_fname = '{:s}_6dist_{:s}_parms.txt'.format(dataset,'comp')
    dist_fpath = os.path.join(BaseDir,dist_fname)
    with open(dist_fpath,'r') as f:
        dist_dict  = json.load(f)
    p_parms_opt = dist_dict['p_parms_opt']
    parms       = dist_dict['parms']
    parms[2] = 'Tau_u'
    
    # load correlations
    DictName   = '{:s}_correlations.mat'.format(dataset)
    DictPath   = os.path.join(BaseDir,DictName)
    CorrDict   = scio.loadmat(DictPath,squeeze_me=True)
    R          = CorrDict['corrs']
    all_fields = [s.rstrip() for s in CorrDict['all_fields']]
    
    # get correlation matrix for that height
    Corr = np.ones((len(ParmSample),len(ParmSample)))
    for iP1 in range(len(ParmSample)):
        parm1 = ParmSample[iP1]
        for iP2 in range(1,len(ParmSample)):
            parm2 = ParmSample[iP2]
            iC1   = iH*len(all_fields) + all_fields.index(parm1)
            iC2   = iH*len(all_fields) + all_fields.index(parm2)
            Corr[iP1,iP2] = R[iC1,iC2]
            Corr[iP2,iP1] = R[iC1,iC2]
    C = np.linalg.cholesky(Corr)
    
    # draw wind samples
    NumSampd  = 0
    iU,iSig,iL,irho = ParmSample.index('Mean_Wind_Speed'),\
        ParmSample.index('Sigma_u'),ParmSample.index('Tau_u'),\
        ParmSample.index('Concentration_u')
    WindParms = np.empty((NumSamps,len(ParmSample)))
    while NumSampd < NumSamps:
        G_unc = np.random.normal(size=(len(ParmSample),NumSamps))
        G_cor = np.dot(C,G_unc)
        U_cor = scipy.stats.norm.cdf(G_cor)
        Y_cor = np.empty(U_cor.shape)
        for iP in range(len(ParmSample)):
            Y_cor[iP,:] = inversecompositeCDF(U_cor[iP,:],
                                                 *p_parms_opt[iP][iH][:-1])
        UInRange = np.logical_and(Y_cor[iU,:] <= URange[1],
                                  Y_cor[iU,:] >= URange[0])
        sigInRange = np.logical_and(Y_cor[iSig,:] <= SigRange[1],
                                    Y_cor[iSig,:] >= SigRange[0])
        LInRange = np.logical_and(Y_cor[iL,:] <= TauRange[1],
                                  Y_cor[iL,:] >= TauRange[0])
        rhoInRange = np.logical_and(Y_cor[irho,:] <= RhoRange[1],
                                    Y_cor[irho,:] >= RhoRange[0])
        InRange = np.logical_and(np.logical_and(np.logical_and(UInRange,
                                                               sigInRange),
                                                LInRange),
                                 rhoInRange)
        Y_cor = Y_cor[:,InRange]
        End   = min(Y_cor.shape[1]+NumSampd,NumSamps)   # end index
        WindParms[NumSampd:End,:] = Y_cor[:,:(End-NumSampd)].T
        NumSampd = End                                  # update number of samps
                                             
    return WindParms


# =============================================================================
# ---------------------------- TURBSIM ANALYSIS -------------------------------
# =============================================================================

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


def WriteTurbSimInputs(InpName,TSDict,TmplDir,WrDir):
    """ Write .inp and .spc (if UserSpec) files given TurbSim parameters
    
        Args:
            InpName (string): name of desired .inp file with extension
            TSDict (dictionary): dictionary with TurbSim parameters
            TmplDir (string): directory with input template files
            WrDir (string): directory to write files to
    """
        
    # template filename
    spctemp = os.path.join(TmplDir,'Template_UsrSpc.spc')
    TStemp  = os.path.join(TmplDir,'Template_TurbSim.inp')
    
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
        fname_spc          = InpName.rstrip('.inp') + '.spc'
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
        
        fpath_spc = os.path.join(WrDir,fname_spc)
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
    fpath_inp = os.path.join(WrDir,InpName)
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
                    f_out.write(line.format('\"'+fpath_spc+'\"'))
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


def MakeTSDict(TurbName,URef,sig_u,L_u,rho_u,R1,R2):
    """ Create TurbSim dictionary for turbine
    """
    import numpy as np
    
    # values that are consistent across all simulations
    T    = 630.0
    dt   = 0.05
    ZRef = 100.0
    mu   = np.pi
    
    # grid parameters 
    ZHub,DZ,DY,n_z,n_y = Turb2Grid(TurbName)
    
    # initialize dictionary
    TS = {}

    #    wind parameters
    TS['ZRef']  = ZRef
    TS['URef']  = URef
    TS['sig_u'] = sig_u
    TS['sig_v'] = 0.8*sig_u
    TS['sig_w'] = 0.5*sig_u
    TS['L_u']   = L_u
    TS['L_v']   = L_u/8.1*2.7
    TS['L_w']   = L_u/8.1*0.66
    TS['rho_u'] = rho_u
    TS['rho_v'] = rho_u
    TS['rho_w'] = rho_u
    TS['mu_u']  = mu
    TS['mu_v']  = mu
    TS['mu_w']  = mu

    #    simulation parameters
    TS['R1']        = R1
    TS['R2']        = R2
    TS['ZHub']      = ZHub
    TS['DY']        = DY
    TS['DZ']        = DZ
    TS['n_z']       = n_z
    TS['n_y']       = n_y
    TS['T']         = T
    TS['dt']        = dt
    TS['TurbModel'] = 'USRINP'
    TS['ProfType']  = 'PL'

    return TS


# =============================================================================
# --------------------------- FAST ANALYSIS -----------------------------------
# =============================================================================

def CreateTurbineDictionary(TurbName,turb_dir,
                            BModes=1,TModes=1,verbose=0):
    """ Convert information in text files to dictionary containing all of the
        turbine parameters necessary to create all FAST files for a WindPACT 
        turbine
        
        Args:
            TurbName (string): turbine name
            turb_dir (string): path to turbine directory
            BModes (Boolean): whether to get blade modes from 
                              Modes v22 output, opt.
            TModes (Boolean): whether to get tower modes from 
                              Modes v22 output, opt.
            
        Returns:
            TurbDict (dictionary): turbine parameters in dictionary
    """
    
    
    if verbose:
        print('\nWriting turbine dictionary' + \
            ' for {:s} from {:s}'.format(TurbName,turb_dir))
    
    # ===================== Initialize turbine dictionary =====================
    TurbDict = {}
    TurbDict['TurbName'] = TurbName
    Date = datetime.datetime.now().strftime('%d-%b-%y')
    TurbDict['FASTCmnt1'] = 'FAST v7.02 input file for turbine \"{:s}\"'.format(TurbName)
    TurbDict['FASTCmnt2'] = 'Generated by J. Rinker (Duke University) ' + \
                        'on {:s} based on WindPACT Excel input files'.format(Date)
                        
    # ======================== Add nacelle properties =========================
    
    # load nacelle data from text file
    fNacName = os.path.join(turb_dir,'parameters\\'+TurbName+'_Nacelle.txt')
    with open(fNacName,'r') as f:
        for line in f:
            row = line.strip('\n').rstrip().split()
            if row:
                key          = row[0]
                value        = float(row[1])
                TurbDict[key] = value
    
    # ========================= Add blade properties ==========================
    
    # load rotor data from text file
    fBldName = os.path.join(turb_dir,'parameters\\'+TurbName+'_Blades.txt')
    with open(fBldName,'r') as f:
        key,i_line = '',0
        while (key != 'DistBldProps'):
            row = f.readline().strip('\n').rstrip().split()
            key = row[0]
            if (key != 'DistBldProps'):
                value = float(row[1])
                for i_bl in range(3):
                    key_bl = '{:s}_{:d}'.format(key,i_bl+1)
                    TurbDict[key_bl] = value
                
            # update line count
            i_line += 1
            
    # calculate blade schedule from remaining table
    n_skip = i_line
    BldSched = np.genfromtxt(fBldName,skip_header=n_skip).tolist()
        
    # load mode shape data from Modes output if available
    fModesName = os.path.join(turb_dir,'modes\\'+TurbName+'_BldModes.mod')
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
    
    # create blade comment
    BldCmnt = 'FAST v7.02 blade file for turbine ' + \
              '\"{:s}\" (J. Rinker, Duke University, {:s})'.format(TurbName,Date)
    
    # save values in dictionary
    for i_bl in range(3):
        TurbDict['BldFile({:d})'.format(i_bl+1)] = TurbName + '_Blade.dat'
        TurbDict['BldSched_{:d}'.format(i_bl+1)] = BldSched
        TurbDict['BldCmnt_{:d}'.format(i_bl+1)]  = BldCmnt
        for i_coeff in range(5):
            TurbDict['BldFl1Sh({:d})_{:d}'.format(i_coeff+2,i_bl+1)] = BldModes[0][i_coeff]
            TurbDict['BldFl2Sh({:d})_{:d}'.format(i_coeff+2,i_bl+1)] = BldModes[1][i_coeff]
            TurbDict['BldEdgSh({:d})_{:d}'.format(i_coeff+2,i_bl+1)] = BldModes[2][i_coeff]
            
    # ========================== Add tower properties =========================
    
    # load tower data from text file
    fTwrName = os.path.join(turb_dir,'parameters\\'+TurbName+'_Tower.txt')
    with open(fTwrName,'r') as f:
        key,i_line = '',0
        while (key != 'TwrSchedFields'):
            row = f.readline().strip('\n').rstrip().split()
            key = row[0]
            if (key != 'TwrSchedFields'):
                value = float(row[1])
                TurbDict[key] = value
                
            i_line += 1
            
    # calculate tower schedule from remaining table
    n_skip                 = i_line
    TwrSched               = np.genfromtxt(fTwrName,skip_header=n_skip).tolist()
    
    # load mode shape data from Modes output if available
    fModesName = os.path.join(turb_dir,'modes',TurbName+'_TwrModes.mod')
    TwrModes = np.empty((2,5))
    TwrModes[:] = np.nan
    if TModes:
        with open(fModesName,'r') as f:
            i_line = 0
            for line in f:
                if ((i_line >= 14) and (i_line <= 18)):
                    for i_mode in range(TwrModes.shape[0]):
                        TwrModes[i_mode,i_line-14] = float(line.split()[i_mode+1])
                i_line += 1
    TwrModes = TwrModes.tolist()
    
    # save variables
    TurbDict['TwrFile']  = TurbName + '_Tower.dat'
    TurbDict['TwrSched'] = TwrSched
    TurbDict['TwrCmnt']  = 'FAST v7.02 tower file for turbine ' + \
                '\"{:s}\" (J. Rinker, Duke University, {:s})'.format(TurbName,Date)
    for i_coeff in range(5):
        TurbDict['TwFAM1Sh({:d})'.format(i_coeff+2)] = TwrModes[0][i_coeff]
        TurbDict['TwFAM2Sh({:d})'.format(i_coeff+2)] = TwrModes[1][i_coeff]
        TurbDict['TwSSM1Sh({:d})'.format(i_coeff+2)] = TwrModes[0][i_coeff]
        TurbDict['TwSSM2Sh({:d})'.format(i_coeff+2)] = TwrModes[1][i_coeff]
           
    # ==================== Add AeroDyn properties ==================
    
    # load AeroDyn data from text file
    fADName = os.path.join(turb_dir,'parameters\\'+TurbName+'_AeroDyn.txt')
    with open(fADName,'r') as f:
        key,i_line = '',0
        while (key != 'ADSchedFields'):
            row = f.readline().strip('\n').rstrip().split()
            key = row[0]
            if (key == 'FoilNm'):
                values = row[1:]
                TurbDict[key] = values
                TurbDict['NumFoil'] = len(values)
            elif (key != 'ADSchedFields'):
                values = float(row[1])
                TurbDict[key] = values
                          
            # save extracted value
            i_line += 1
            
        # calculate AeroDyn schedule from remaining table
        ADSched = []
        for i_AD in range(int(TurbDict['BldNodes'])):
            row = f.readline().strip('\n').rstrip().split()
            ADSchedRow = [float(x) for x in row[:-1]] + [row[-1]]
            ADSched.append(ADSchedRow)
    
    # manually calculate hub height for AD
    TurbDict['HH']      = TurbDict['TowerHt'] + TurbDict['Twr2Shft'] + \
                            TurbDict['OverHang']*np.sin(TurbDict['ShftTilt']*np.pi/180)
    
    # save variables
    TurbDict['ADFile']  = TurbName + '_AD.ipt'
    TurbDict['ADSched'] = ADSched
    TurbDict['ADCmnt']  = 'FAST v7.02 AeroDyn file for turbine ' + \
                '\"{:s}\" (J. Rinker, Duke University, {:s})'.format(TurbName,Date)
                
    # ======================== Add control properties =========================
    
    # load control data from text file
    fCntrName = os.path.join(turb_dir,'parameters\\'+TurbName+'_Control.txt')
    with open(fCntrName,'r') as f:
        for line in f:
            row = line.strip('\n').rstrip().split()
            if row:
                key = row[0]
                if ('P2P' in key):
                    value = [float(s) for s in row[1:]]
                else:
                    value = float(row[1])
                TurbDict[key] = value
                
    # save variables
    TurbDict['PitchCmnt']  = 'Pitch controller for turbine \"{:s}\"'.format(TurbName) + \
                                ' (J. Rinker, Duke University, {:s})'.format(Date)
                
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


def writeBldModes(fpath_temp,fpath_out,TurbDict,BldIdx):
    """ Blade input file for Modes v22
    """
    
    
    # get interpolated blade structural parameters
    BldSched = TurbDict['BldSched_{:.0f}'.format(BldIdx)]
    RatedRPM = TurbDict['CNST(2)']
    TipRad   = TurbDict['TipRad']
    HubRad   = TurbDict['HubRad']
        
    # open template file and file to write to
    with open(fpath_temp,'r') as f_temp:
        with open(fpath_out,'w') as f_write:
            i_line = 0
            for line in f_temp:
                if i_line == 1:
                    f_write.write(line.format(RatedRPM))
                elif i_line == 3:
                    f_write.write(line.format(TipRad))
                elif i_line == 4:
                    f_write.write(line.format(HubRad))
                elif i_line == 5:
                    f_write.write(line.format(0.0))
                elif i_line == 8:
                    f_write.write(line.format(len(BldSched)))
                elif i_line == 12:
                    for i_BlNode in range(len(BldSched)):
                        row = [BldSched[i_BlNode][0],0.,
                               BldSched[i_BlNode][3],
                               BldSched[i_BlNode][4],
                               BldSched[i_BlNode][5]]
                        f_write.write(line.format(*row))
                else:
                    f_write.write(line)
                i_line += 1
    
    return
    

def writeTwrModes(fpath_temp,fpath_out,TurbDict):
    """ Tower input file for Modes v22
    """
    
    # load tower-modes-specific vales
    TwrSched   = TurbDict['TwrSched']
    TowerHt       = TurbDict['TowerHt']
    TTMass        = TurbDict['NacMass'] + TurbDict['HubMass']
    
    # open template file and file to write to
    with open(fpath_temp,'r') as f_temp:
        with open(fpath_out,'w') as f_write:
            i_line = 0
            for line in f_temp:
                if i_line == 3:
                    f_write.write(line.format(TowerHt))
                elif i_line == 5:
                    f_write.write(line.format(TTMass))
                elif i_line == 8:
                    f_write.write(line.format(len(TwrSched)))
                elif i_line == 12:
                    for i_TwrNode in range(len(TwrSched)):
                        row = [TwrSched[i_TwrNode][0],
                               TwrSched[i_TwrNode][1],
                               TwrSched[i_TwrNode][2],
                               TwrSched[i_TwrNode][3]]
                        f_write.write(line.format(*row))
                else:
                    f_write.write(line)
                i_line += 1
    
    return

def writeBlade(fpath_temp,fpath_out,TurbDict):
    """ Blade input file for FAST v7.02
    """
    
    
    # calculate blade-specific vales
    TurbName = TurbDict['TurbName']
    Comment  = 'FAST v7.02 blade file for turbine \"{:s}\"'.format(TurbName) + \
                ' (JRinker, Duke University)'
    Damping  = TurbDict['BldDmp']
    BldSched = TurbDict['BldSched']
    BldModes = np.array(TurbDict['BldModes'])
    
    # open template file and file to write to
    with open(fpath_temp,'r') as f_temp:
        with open(fpath_out,'w') as f_write:
            i_line = 0
            for line in f_temp:
                if i_line == 2:
                    f_write.write(line.format(Comment))
                elif i_line == 4:
                    f_write.write(line.format(len(BldSched)))
                elif i_line == 6:
                    f_write.write(line.format(Damping[0]))
                elif i_line == 7:
                    f_write.write(line.format(Damping[1]))
                elif i_line == 8:
                    f_write.write(line.format(Damping[2]))
                elif i_line == 18:
                    for i_BlNode in range(len(BldSched)):
                        f_write.write(line.format(*BldSched[i_BlNode]))
                elif ((i_line >= 20) and (i_line <= 34)):
                    f_write.write(line.format(
                        BldModes.reshape(BldModes.size)[i_line-20]))
                else:
                    f_write.write(line)
                i_line += 1
    
    return
    
    
def writeTower(fpath_temp,fpath_out,TurbDict):
    """ Tower input file for FAST v7.02
    """
    
    
    # calculate tower-specific vales
    TurbName = TurbDict['TurbName']
    Comment  = 'FAST v7.02 tower file for turbine \"{:s}\"'.format(TurbName) + \
                ' (JRinker, Duke University)'
    Damping  = TurbDict['TwrDmp']
    TwrModes = np.array(TurbDict['TwrModes'])
    TwrSched = np.array(TurbDict['TwrSched'])
        
    # open template file and file to write to
    with open(fpath_temp,'r') as f_temp:
        with open(fpath_out,'w') as f_write:
            i_line = 0
            for line in f_temp:
                if i_line == 2:
                    f_write.write(line.format(Comment))
                elif i_line == 4:
                    f_write.write(line.format(len(TwrSched)))
                elif i_line == 6:
                    f_write.write(line.format(Damping[0]))
                elif i_line == 7:
                    f_write.write(line.format(Damping[1]))
                elif i_line == 8:
                    f_write.write(line.format(Damping[2]))
                elif i_line == 9:
                    f_write.write(line.format(Damping[3]))
                elif i_line == 21:
                    for i_TwrNode in range(len(TwrSched)):
                        f_write.write(line.format(*TwrSched[i_TwrNode,:]))
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
    
    
def WriteAeroDynTemplate(fpath_temp,fpath_out,TurbDict):
    """ AeroDyn input file for FAST v7.02
    """
    
    # calculate file-specific vales
    TurbName       = TurbDict['TurbName']
    Comment1       = 'AeroDyn v13.00 file for turbine ' + \
                        '\"{:s}\" (JRinker, Duke University)'.format(TurbName)
    WindHH         = TurbDict['TowerHt'] + TurbDict['Twr2Shft'] + \
                        TurbDict['OverHang']*np.sin(TurbDict['ShftTilt']*np.pi/180.)
    NumFoil        = len(TurbDict['FoilNm'])
    
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
                    
                    # special checks for some fields
                    if (field == 'FoilNm'):
                        AF_fname = TurbDict[field][0]+'.dat'
                        AF_fpath = os.path.join('AeroData',AF_fname)
                        w_line = r_line.format(AF_fpath)
                        f_write.write(w_line)
                        for i_foil in range(1,NumFoil):
                            AF_fname = TurbDict[field][i_foil]+'.dat'
                            AF_fpath = os.path.join('AeroData',AF_fname)
                            w_line = '\"{:s}\"\n'.format(AF_fpath)
                            f_write.write(w_line)
                            
                    elif (field == 'BldNodes'):
                        w_line = r_line.format(TurbDict[field])
                        f_write.write(w_line)
                        f_write.write('RNodes    AeroTwst  ' + \
                                        'DRNodes  Chord  NFoil  PrnElm\n')
                        for i_ADnode in range(int(TurbDict[field])):
                            w_line = '{:8.5f}{:7.2f}{:12.5f}{:7.3f}{:3.0f}'.format( \
                                        *TurbDict['ADSched'][i_ADnode]) \
                                        + '      NOPRINT\n'
                            f_write.write(w_line)
                        return
                            
                    else:
                    
                        # try to load value from turbine dictionary
                        try:
                            w_line = r_line.format(TurbDict[field])
                            
                        # if key isn't present, check a few conditions
                        except KeyError:
                            if (field == 'HH'):
                                w_line = r_line.format(WindHH)
                            elif (field == 'NumFoil'):
                                w_line = r_line.format(NumFoil)
                            else:
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
                        
                        # check subkey
                        subkey = [sk for sk in DictFields if sk in field] 
                        if subkey:
                            w_line = r_line.format(TurbDict[subkey[0]])
                            
                        # check if it's an unused file
                        elif (('File' in field) and (field != 'ADFile')):
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
                   wind_dir=None,fileID=None,
                   wr_dir=None,
                   **kwargs):
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
    
    
    
    
    # get initial wind speed
    wind_fpath = os.path.join(wind_dir,wind_fname)
    u0 = GetFirstWind(wind_fpath)
    
    print('Writing FAST files for \"{:s}\" '.format(TName) + \
            'with wind file {:s}'.format(wind_fpath))
            
    # default ICDict
    ICDict = {}
    ICDict['BlPitch(1)'],ICDict['BlPitch(2)'],ICDict['BlPitch(3)'] = 2.6,2.6,2.6
    ICDict['OoPDefl'],ICDict['IPDefl'],ICDict['TeetDefl'] = 0.,0.,0.
    ICDict['Azimuth'],ICDict['RotSpeed'],ICDict['NacYaw'] = 0.,6.,0.
    ICDict['TTDspFA'],ICDict['TTDspSS'] = 0.,0.
    ICDict['GenDOF'],ICDict['TMax'],ICDict['TStart'] = 'True',630.,30.
    
    # set optional values as necessary
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
        
    # load passed-in ICs
    for key in kwargs.keys():
        if key in ICDict.keys():
            ICDict[key] = kwargs[key]

    # create filenames
    fAD_temp  = os.path.join(turb_dir,'templates',TName+'_AD_template.ipt')
    fAD_out   = os.path.join(wr_dir,fAD_name)
    fFST_temp = os.path.join(turb_dir,'templates',TName+'_template.fst')
    fFST_out  = os.path.join(wr_dir,fFST_name)
    
    ICDict['ADFile'] = fAD_out
    
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
                if ('{:' in line):
                    key    = line.split()[1]
                    valfmt = line.split()[0]
                    cmnt  = line.split(valfmt)[1]
                    if key in ICDict.keys():
                        valstr = valfmt.format(ICDict[key])
                        line = ''.join([valstr,cmnt])
                    else:
                        print('**** NO VALUE FOR KEY {:s} ****'.format(key))
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
    
    # define "global" parameters for simplicity
    legFS = 'x-small'
    
    # create figure if none specified
    if fig is None:
        fig = plt.figure(figsize=(7,9))
        
    # axes locations
    xPlot = 0.11
    yPlot = np.arange(3)[::-1]*0.32 + 0.10
    wd,ht = 0.80,0.18
    
    # Axes 1: GenSpeed,RotPwr,GenPwr,RotThrust,RotTorq
    ax1 = fig.add_axes([xPlot,yPlot[0],wd,ht])
    
    plt.plot(x,Data[:,Fields.index('GenSpeed')],label='GenSpeed (rpm)')
#    plt.plot(x,Data[:,Fields.index('LSShftPwr')],label='RotPwr (kW)')
    plt.plot(x,Data[:,Fields.index('GenPwr')],label='GenPwr (kW)')
    plt.plot(x,Data[:,Fields.index('RotThrust')],label='RotThrust (kN)')
    plt.plot(x,Data[:,Fields.index('RotTorq')],label='RotTorq (kN-m)')
    
    ax1.set_xlim([x[0],x[-1]])
    ax1.grid('on')
#    ax1.legend(loc=2,fontsize=legFS)
    ax1.legend(fontsize=legFS,
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode='expand', borderaxespad=0.)
    ax1.set_xticks(np.arange(x[0],x[-1]+1,2))
    
    
    # Axes 2: RotSpeed,BlPitch,GenTq,TSR
    ax2 = fig.add_axes([xPlot,yPlot[1],wd,ht])
    
    plt.plot(x,Data[:,Fields.index('RotSpeed')],label='RotSpeed (rpm)')
    plt.plot(x,Data[:,Fields.index('BldPitch1')],label='BlPitch ($^\mathrm{o}$)')
    plt.plot(x,Data[:,Fields.index('GenTq')],label='GenTq (kN-m)')
    plt.plot(x,Data[:,Fields.index('TSR')],label='TSR (-)')
    
    ax2.set_xlim([x[0],x[-1]])
    ax2.grid('on')
#    ax2.legend(loc=2,fontsize=legFS)
    ax2.legend(fontsize=legFS,
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode='expand', borderaxespad=0.)
    ax2.set_xticks(np.arange(x[0],x[-1]+1,2))
    
    # Axes 3: OopDefl1,IPDefl1,TTDspFA,TTDspSS
    ax3 = fig.add_axes([xPlot,yPlot[2],wd,ht])
    
    plt.plot(x,Data[:,Fields.index('OoPDefl1')],label='OoPDefl1 (m)')
    plt.plot(x,Data[:,Fields.index('IPDefl1')],label='IPDefl1 (m)')
    plt.plot(x,Data[:,Fields.index('TTDspFA')],label='TTDspFA (m)')
    plt.plot(x,Data[:,Fields.index('TTDspSS')],label='TTDspSS (m)')
    
    #ax3.set_ylim([0,40])
    ax3.set_xlim([x[0],x[-1]])
    ax3.grid('on')
#    ax3.legend(loc=2,fontsize=legFS)
    ax3.legend(fontsize=legFS,
               bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode='expand', borderaxespad=0.)
    ax3.set_xticks(np.arange(x[0],x[-1]+1,2))
    ax3.set_xlabel('Wind speed (m/s)',
                   labelpad=10)
    
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
    
    
def CalculateDELs(x,m,
                 lookahead=2,EqFreq=1.,T_DEL=600.,PSF=1.):
    """ Damage equivalent load for time series with given SN slopes
    
        Args:
            x (numpy.ndarray): time series
            m (float/numpy.ndarray/list): SN-slopes
            lookahead (int): steps to look ahead in peak detect [opt]
            EqFreq (float): equivalent frequency, Hz [opt]
            T_DEL (float): time span of DEL calc, s [opt]
            PSF (float): partial safety factor [opt]
            
        Returns:
            DEL (float/numpy.ndarray): damage equivalent loads
    """
    
    # change m to 1D array if float passed in
    InFloat = 0
    if (type(m) in ['float']):
        m = np.array(m)
        InFloat = 1
    
    # get extremal values from time series
    maxinfo, mininfo = peakdetect(x, 
                                  lookahead=lookahead)
    ipeaks = np.array([tup[0] for tup in maxinfo] + [tup[0] for tup in mininfo])
    xpeaks = np.array([tup[1] for tup in maxinfo] + [tup[1] for tup in mininfo])
    
    # re-sort algorithm into increasing indices
    isort  = ipeaks.argsort()
    ipeaks = ipeaks[isort]
    xpeaks = xpeaks[isort]
    
    # do rainflow counting
    rf_counts = rainflow(xpeaks)
    
    # calculate DELs
    N_eq = T_DEL * EqFreq                       # equivalent cycles
    n, S = rf_counts[3,:], rf_counts[2,:]       # cycle counts and amplitudes
    DEL = np.empty(len(m))
    for i_m in range(len(m)):                   # loop through slope values
        DEL[i_m] = ( np.sum(n * ( S )**m[i_m]) / ( N_eq ) ) ** \
                        (1. / m[i_m]) * PSF
    
    # turn back to float if that's what was given
    if InFloat: DEL = float(DEL)
        
    return DEL
    

def CalculateDELsAll(FastPath,
                     m_comp=[8.,10.,12.],m_metal=[3.,4.,5.],
                     t_lookahead=0.1):
    """ Calculate all DELs for FAST file
    
        Args:
            FastPath (string): path to output FAST file
            m_comp (float/list/np.ndarray): SN slopes for composite material [opt]
            m_metal (float/list/np.ndarray): SN slopes for metal material [opt]
            t_lookahead (float): lookahead time for peak detect algorithm [opt]
            
        Returns:
            DELDict (dictionary): dictionary with DEL information
    """
    
    # load FAST file
    FastDict = ReadFASTFile(FastPath)
    fields = FastDict['Fields']
    units  = FastDict['Units']
    data   = FastDict['Data']
    
    # get lookahead in number of steps
    time = data[:,fields.index('Time')]
    dt = (time[2] - time[0]) / 2.
    lookahead = int(np.round(t_lookahead / (dt)))
    
    # get list of keys that are forces/moments
    idcs_FM = [i for i in range(len(fields)) if ('N' in units[i])]
    fields_FM = [fields[i] for i in idcs_FM]
    
    # iteratively calculate DELs
    DELkeys  = []
    DELunits = []
    DELs     = []
    for i_key in range(len(idcs_FM)):
        
        # get key, time series, and key units
        key    = fields_FM[i_key]
        x      = data[:,fields.index(key)]
        kunits = units[fields.index(key)]
        
        # key corresponds to composite item
        if any([s in key for s in ('Root','Spn')]):
            
            # get DELs for all slopes
            m_all = np.array(m_comp)
            DELs_key = CalculateDELs(x,m_all,
                                      lookahead=lookahead)
                                      
            # save DEL key, units, and value
            for i_m in range(len(m_all)):
                m = m_all[i_m]
                DELkey = '{:s}_m={:.1f}'.format(key,m)
                DELkeys.append(DELkey)
                DELunits.append(kunits)
                DELs.append(DELs_key[i_m])
                                
        # key corresponds to metal item
        elif any([s in key for s in ('LSS','HSS','YawBr','Twr')]):
            
            # get DELs for all slopes
            m_all = np.array(m_metal)
            DELs_key = CalculateDELs(x,m_all,
                                      lookahead=lookahead)
                                      
            # save DEL key, units, and value
            for i_m in range(len(m_all)):
                m = m_all[i_m]
                DELkey = '{:s}_m={:.1f}'.format(key,m)
                DELkeys.append(DELkey)
                DELunits.append(kunits)
                DELs.append(DELs_key[i_m])
                
    # save DEL data in dictionary
    DELDict = {}
    DELDict['Keys']  = DELkeys
    DELDict['Units'] = DELunits
    DELDict['DELs']   = DELs
    
    return DELDict
                         

def polyregression(x,y,p_i):
    """ Fit multi-dimensional polynomial surface to input data and output data
        given individual polynomial powers p_i and excluding cross terms with
        combined power larger than max(p_i).
        
        Args:
            x (list/numpy array): input data
            y (list/numpy array): output data
            p_i (list/numpy array): highest powers of surface for input vars
            
        Returns:
            coeffs (numpy array): 
    """
    import numpy as np
    
    # convert all input to numpy arrays in case list is fed in
    x   = np.array(x)
    y   = np.array(y)
    p_i = np.array(p_i)
    
    # if only 1D in x, reshape to 2D column array
    if len(x.shape) == 1:
        x = x.reshape((x.size,1))
                
    # get max power and check dimension compatibility
    n_x = x.shape[1]
    n_y = x.shape[0]    
    if (y.size != n_y):
        errStr = 'Dimensions of x and y are not compatible'
        raise ValueError(errStr)
    if (len(p_i) != n_x):
        errStr = 'Dimensions of x and p_i are not compatible'
        raise ValueError(errStr)
        
    # get list of powers
    ps = GetAllPowers(p_i)
          
    # create Vandermonde matrix
    A = myvander(x,ps)
        
    # solve least squares problem to get coefficients
    coeffs = np.linalg.lstsq(A,y)[0]
    
    return (coeffs,ps)
    

def OLSfit(Xv, y):
    """ OLS fit to data in Xv and y
    
        Args:
            x (numpy array): array of input data
            y (numpy array): array of output data
            
        Returns:
            results (RegressionResultsWrapper): output from sm.OLS.fit
    """
    
    # solve OLS problem
    results = sm.OLS(y, Xv).fit()
    
    return results  
    

def SigOLSfit(x, y, pmax_i,
              alpha=0.05):
    """ Significant OLS fit to data in x and y with max poly order pmax_i
    
        Args:
            x (numpy array): array of input data
            y (numpy array): array of output data
            p_i (list/numpy array): list of polynomial orders
            alpha (float): significance threshold
            
        Returns:
            results (RegressionResultsWrapper): output from sm.OLS.fit
    """
    
    # solve basic OLS problem
    results = OLSfit(x, y)
    
    # get significant coefficients
    ps_all = GetAllPowers(pmax_i)
    ps_red = ps_all[results.pvalues <= alpha]
    X_red  = x[:,results.pvalues <= alpha]
    cs_red = OLSfit(X_red, y).params
    
    return X_red, cs_red, ps_red    
    

def GetAllPowers(p_i):
    """ Array of all powers for maximum powers in p_i
    
        Args:
            p_i (list/numpy array): maximal power for each individual term
            
        Returns:
            ps (list/numpy array): list of powers with cross-terms
    """
    
    # convert input to numpy arrays in case list is fed in
    p_i = np.array(p_i)
           
    # get max power
    p_max = p_i.max()

    # get all possible lists of powers with cross terms
    n_all = np.prod(p_i+1.)                 # total number of poss. power combos
    p = []                                  # list of ind. powers (0,1,2,...)
    for i in range(len(p_i)):
        p.append(np.arange(p_i[i]+1))
    P = np.meshgrid(*p)                     # arrays of cross-powers
    ps = np.empty((n_all,len(p_i)))         # convert arrays to vectors
    for i in range(len(p_i)):
        ps[:,i] = P[i].reshape(n_all)       # save vectors in total array

    # remove combinations of powers with total power > p_max
    ps = ps[np.sum(ps,axis=1) <= p_max] 
        
    return ps
    
    
def myvander(x,ps):
    """ A sort-of Vandermonde matrix with the powres defined in array ps. Used
        in polyregression routine.
        
        Args:
            x (numpy array): input data
            ps (numpy array): array of powers
            
        Returns:
            A (numpy array): 2D array of x raised to powers in ps
    """
    import numpy as np
    
    # if only 1D in x or p, reshape to 2D array
    if len(x.shape) == 1:
        x = x.reshape((x.size,1))
    if len(ps.shape) == 1:
        ps = ps.reshape((1,ps.size))
        
    # define dimension parameters
    n_y = x.shape[0]
    n_coeff = ps.shape[0]        
        
    # creaate matrix
    A = np.empty((n_y,n_coeff))
    for i in range(n_coeff):
        A[:,i] = np.prod(x ** ps[i,:], axis=1)
        
    return A
    
    

def DiscreteOpt(ErrFnc,p0,
                max_iters=1000,err_thresh=1e-10,
                derr_thresh=0.01,n_last=4,
                verbose=0):
    """ Optmize discrete error given initial guess p0
    
        Args:
            ErrFnc (function): function to return error given guess p
            p0 (list/np.array): initial guess
            max_iters (int): maximum number of iterations [opt]
            err_thresh (float): minimum threshold for error [opt]
            derr_thresh (float): threshold for error diff over last few iters [opt]
            n_last (int): number of iterations to look back, see if decreasing [opt]
            verbose (int): flag to print iteration updates
            
        Returns:
            results (dictionary): optimization results
    """
    
    print('\nRunning optimization...')
    
    # initialize variables
    p          = p0                     # array of optimization variables
    iter_count = 1                      # iteration counter
    stop       = 0                      # flag to stop iterations
    curr_err   = ErrFnc(p)              # current error value
    last_errs  = (1e16)*np.ones(n_last) # errors for last several iterations
    n_x        = len(p0)                # number of optimization variables
    last_ps    = np.empty((n_last,n_x)) # last set of guesses
    results    = {}                     # dictionary for output
    
    # iterate until converge
    last_ps[0,:] = p0
    while(not stop):
        
        if verbose:
            print('Iteration {:d}. Error = {:.2f}'.format(iter_count,curr_err) + \
                    ' p_init = [{:s}]'.format(','.join([str(x) for x in p])))
            
        # iteratively increment each variable, calculating and saving error    
        iter_errs = np.empty(n_x)
        for i_x in range(n_x):
            ptest = np.copy(p)
            ptest[i_x] += 1
            iter_errs[i_x] = ErrFnc(ptest)
            
        # determine variable to increment
        i_incr = iter_errs.argmin()
        min_err = iter_errs[i_incr]
        
        if verbose:
            print('  incrementing x{:d} (new_err = {:.1f})'.format(i_incr,min_err))
                                 
        # if new min error is 5% larger than current error, save current and stop
        if (min_err >= 1.05*curr_err):
            print('\nOptimization halting: local minimum achieved\n')
            results['num_iters'] = iter_count
            results['p_out']     = p
            results['error']     = curr_err
            results['exit_code'] = 0
            return results
            
        # if new min error is smaller than error threshold
        if (min_err <= err_thresh):
            print('\nOptimization halting: error is below threshold\n')
            results['num_iters'] = iter_count
            results['p_out']     = p
            results['error']     = curr_err
            results['exit_code'] = 0
            return results
            
        # if there is only a small change in error for the last few iterations
        if (np.abs((curr_err - last_errs[-1])/curr_err) <= derr_thresh):
            print('\nOptimization halting: no large reduction in error\n')
            results['num_iters'] = iter_count - n_last
            results['p_out']     = last_ps[-1,:]
            results['error']     = last_errs[-1]
            results['exit_code'] = 0
            return results
            
        # if max number of iterations reached
        if iter_count >= max_iters:
            print('\nOptimization halting: maximum iteration count reached\n')
            results['num_iters'] = iter_count
            results['p_out']     = p
            results['error']     = curr_err
            results['exit_code'] = 1
            return results
            
        # update values
        p[i_incr]    += 1
        curr_err      = min_err
        iter_count   += 1
        last_errs[1:] = last_errs[:-1]
        last_errs[0]  = curr_err
        last_ps[1:,:] = last_ps[:-1,:]
        last_ps[0,:]  = p
            
    return results
    
    
def FASTUniqueStats(x,y,RunName):
    """ Calculate the mean and standard deviation across duplicates in vector
        y as determined by the input data in x
        
        Args:
            x (numpy array): input data
            y (numpy array): output data
            RunName (string): which run to analyze
            
        Returns:
            x_uniq (numpy array): unique input points
            err_mean (numpy array): mean of error
            err_std (numpy array): std dev of error
    """
    
    # get list of wind parameters for that run
    WindParmDict = RunName2WindParms(RunName)
    URefs, Is, Ls, rhos = WindParmDict['URefs'],WindParmDict['Is'], \
                          WindParmDict['Ls'],WindParmDict['rhos']
    WindParms = [URefs,Is,np.log10(Ls),rhos]
    
    # loop through wind parameters
    n_uniq = np.prod(np.array([len(l) for l in WindParms]))
    x_uniq = np.empty((n_uniq,x.shape[1]))
    y_mean = np.empty(n_uniq)
    y_std  = np.empty(n_uniq)
    i_uniq   = 0
    for iU,iI,iL,iRho in [(a,b,c,d) for a in range(len(WindParms[0])) \
                                    for b in range(len(WindParms[1])) \
                                    for c in range(len(WindParms[2])) \
                                    for d in range(len(WindParms[3]))]:
        
        U,I,logL,rho = WindParms[0][iU], WindParms[1][iI], \
                      WindParms[2][iL], WindParms[3][iRho]
                      
        mask = np.logical_and(np.logical_and(np.logical_and(x[:,0]==U,
                                                            x[:,1]==I),
                                             x[:,2]==logL),
                              x[:,3]==rho)

        # assign data to aray
        x_uniq[i_uniq] = [U,I,logL,rho]
        y_data         = y[mask]
        y_mean[i_uniq] = np.mean(y_data)
        y_std[i_uniq]  = np.std(y_data)
        
        i_uniq += 1
        
    return x_uniq, y_mean, y_std

# =============================================================================
# ---------------------------------- MAPPINGS ---------------------------------
# =============================================================================
    
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
    
    if (dataset in ('NREL','fluela')):

        # extract date information
        year     = int(fname[6:10])
        month    = int(fname[0:2])
        day      = int(fname[3:5])
        hour     = int(fname[11:13])
        minute   = int(fname[14:16])
        time_tup = (year,month,day,hour,minute)
    
        # convert tuple to float
        time_flt = timetup2flt(time_tup)
        
    elif (dataset == 'PM06'):

        # extract date information
        year     = int(fname[6:10])
        month    = int(fname[0:2])
        day      = int(fname[3:5])
        hour     = int(fname[11:13])
        minute   = int(fname[13:15])
        time_tup = (year,month,day,hour,minute)
    
        # convert tuple to float
        time_flt = timetup2flt(time_tup)
        
    elif (dataset == 'texastech'):

        # extract date information
        date_str = re.search('(D[0-9]{8})',fname).group()[1:]
        time_str = re.search('(T[0-9]{4})',fname).group()[1:]
        year     = int(date_str[:4])
        month    = int(date_str[4:6])
        day      = int(date_str[6:])
        hour     = int(time_str[:2])
        minute   = int(time_str[2:4])
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
    
    # get tuple of timestampe
    if isinstance(timestamp,tuple):
            time_tup = timestamp
    elif isinstance(timestamp,(int,float)):
            time_tup = timeflt2tup(timestamp)
    else:
        errStr = 'Invalid {} for timestamp'. \
                     format(type(timestamp))
        raise TypeError(errStr)

    if (dataset in ['NREL','fluela']):

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

    elif (dataset== 'PM06'):

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
        fname = '_'.join([monthS,dayS,yearS,hourS+minS,'TS','WND'])
        fpath = os.path.join(dirpath,fname)+'.mat'

        # check if file doesn't exist
        if len(fpath) > 0:
            fpath = fpath[0]
        else:
            errStr = 'File {} does not exist.'.format(fname_part)
            raise IOError(errStr)

    elif (dataset== 'texastech'):

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
        fname_part = '*D'+yearS+monthS+dayS+'_T'+hourS+minS
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
                             
    elif dataset == 'PM06':
        
        grp_nmbr,inst_nmbr = id2grpnmbr(dataset,ID)

        if   (field == 'Sonic_x'):     datfield = 'Sonic_x_{:d}'.format(ID)
        elif (field == 'Sonic_y'):     datfield = 'Sonic_y_{:d}'.format(ID)
        elif (field == 'Sonic_z'):     datfield = 'Sonic_z_{:d}'.format(ID)
        elif (field == 'Sonic_T'):     datfield = 'Sonic_T_{:d}'.format(ID)
        elif (field == 'Humidity'):    datfield = 'Hygro_{:d}'.format(ID)
        elif (field == 'Grp_Warning'): datfield = 'GrpWarn_{:d}'.format(grp_nmbr)
        else:
            errStr = 'Unknown custom field {} for dataset \"{}\"'.format(field,dataset)
            raise AttributeError(errStr)
            
# TODO: uncomment check datfield
            
        # check that datafield is valid
#        if (not check_datfields(dataset,datfield)):
#            raise ValueError('Invalid identifier {:d} for field {:s}'\
#                             .format(ID,field))

    elif dataset == 'texastech':

        if   (field == 'UVW_x'):          datfield = 'UVW_x_{:.0f}m'.format(ID)
        elif (field == 'UVW_y'):          datfield = 'UVW_y_{:.0f}m'.format(ID)
        elif (field == 'UVW_z'):          datfield = 'UVW_z_{:.0f}m'.format(ID)
        elif (field == 'Sonic_x'):        datfield = 'Sonic_x_{:.0f}m'.format(ID)
        elif (field == 'Sonic_y'):        datfield = 'Sonic_y_{:.0f}m'.format(ID)
        elif (field == 'Sonic_z'):        datfield = 'Sonic_z_{:.0f}m'.format(ID)
        elif (field == 'Sonic_T'):        datfield = 'Sonic_T_{:.0f}m'.format(ID)
        elif (field == 'Temperature'):    datfield = 'Air_Temp_{:.0f}m'.format(ID)
        elif (field == 'Humidity'):       datfield = 'Rel_Hum_{:.0f}m'.format(ID)
        elif (field == 'Pressure'):       datfield = 'Baro_Presr_{:.0f}m'.format(ID)
        else:
            errStr = 'Unknown custom field {} for dataset \"{}\"'.format(field,dataset)
            raise AttributeError(errStr)
            
        # check that datafield is valid
        if (not check_datfields(dataset,datfield)):
            raise ValueError('Invalid height {} for field {}'.format(\
                ID,field))

    else:
        errStr = 'Dataset \"{}\" not coded.'.format(dataset)
        raise KeyError(errStr)

    return datfield


def rawfield2datfield(dataset,rawfield):
    """ Convert field from raw data to field used in high-frequency .mat

        Args:
            dataset (string): flag for dataset
            rawfield (string): raw data field

        Returns:
            datfield (string): dataset-specific fieldname for high-frequency data
    """


    if dataset == 'PM06':

        if ('U' in rawfield):
            ID = 4*(int(rawfield[3])-1) + int(rawfield[5])
            datfield = 'Sonic_{:s}_{:.0f}'.format(rawfield[1],ID)
            
        elif ('Ts' in rawfield):
            ID = 4*(int(rawfield[3])-1) + int(rawfield[5])
            datfield = 'Sonic_T_{:.0f}'.format(ID)
            
        elif ('kh2o' in rawfield):
            ID = rawfield[5]
            datfield = 'Hygro_{:s}'.format(ID)
            
        elif ('grp' in rawfield):
            ID = rawfield[-4]
            datfield = 'GrpWarn_{:s}'.format(ID)
            
        elif ('TIME' in rawfield):
            datfield = 'Time'

        else:
            errStr = 'Unknown raw fieldname {:s} '.format(rawfield) + \
                    'for dataset \"{:s}\"'.format(dataset)
            raise AttributeError(errStr)
            
#        # check that datafield is valid
#        if (not check_datfields(dataset,datfield)):
#            raise ValueError('Invalid height {} for field {}'.format(\
#                ID,field))

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
    
    if (dataset == 'PM06'):
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


def RunName2WindParms(RunName):
    """ Simulation parameters for given run name
        
        Args:
            RunName (string): key for run of interest
            
        Returns:
            Parms (dictionary): dictionary of wind parameters
    """
    
    Parms = {}    
    if (RunName in ('Peregrine','BigRun2')):
        Parms['URefs']  = [5,7,9,10,10.5,11,11.5,12,13,16,19,22]
        Parms['Is']     = [0.1,0.2,0.3,0.4,0.5]
        Parms['Ls']     = [10**1.5,10**2.,10**2.5,10**3]
        Parms['rhos']   = [0.,0.1,0.2,0.3,0.4]
        Parms['n_dups'] = 5
        Parms['TurbNames'] = ['WP0.75A08V00','WP1.5A08V03', \
                              'WP3.0A02V02','WP5.0A04V00']
    elif (RunName == 'TestRun'):
        Parms['URefs']  = [7.0, 9.0, 13.0]
        Parms['Is']     = [0.3]
        Parms['Ls']     = [10**2.0]
        Parms['rhos']   = [0.4]
        Parms['n_dups'] = 10
    elif (RunName == 'SmallRun'):
        Parms['URefs']  = [7.0, 9.0, 10.0, 11.0, 13.0, 19.0]
        Parms['Is']     = [0.1,0.3,0.5]
        Parms['Ls']     = [10**2.0]
        Parms['rhos']   = [0.4]
        Parms['n_dups'] = 10
        Parms['TurbNames'] = ['WP0.75A08V00','WP1.5A08V03', \
                              'WP3.0A02V02','WP5.0A04V00']
    elif (RunName == 'TestBig'):
        Parms['URefs']  = [5,7,9,10]
        Parms['Is']     = [0.1,0.2]
        Parms['Ls']     = [10**1.5]
        Parms['rhos']   = [0.]
        Parms['n_dups'] = 1
        Parms['TurbNames'] = ['WP5.0A04V00']
    elif (RunName == 'Fine'):
        Parms['URefs']  = [5,15]
        Parms['Is']     = [0.2]
        Parms['Ls']     = [10**2.5]
        Parms['rhos']   = [0.,0.1,0.3]
        Parms['n_dups'] = 1000
        Parms['TurbNames'] = ['WP5.0A04V00']
        
    return Parms

def FastParm2Name(Parm):
    """ Convert FAST key to plotting name and units
    
        Args:
            Parm (string): FAST parameter
            
        Returns:
            Name (string): plotting name
            Units (string): units
    """
    
    if Parm == 'WindVxi':       Name,Units = 'Hub-Height u','m/s'
    elif Parm == 'WindVyi':     Name,Units = 'Hub-Height v','m/s'
    elif Parm == 'WindVzi':     Name,Units = 'Hub-Height w','m/s'
    elif Parm == 'RotSpeed':    Name,Units = 'Rotor S[eed','rpm'
    elif Parm == 'GenTq':       Name,Units = 'Generator Torque','kN-m'
    elif Parm == 'GenPwr':      Name,Units = 'Generated Power','kW'
    elif Parm == 'TwrBsMyt':    Name,Units = 'Tower Base Fore-Aft Moment','kN-m'
    elif Parm == 'RootMIP1':    Name,Units = 'Bld 1 IP Root Moment','kN-m'
    elif Parm == 'RootMIP2':    Name,Units = 'Bld 2 IP Root Moment','kN-m'
    elif Parm == 'RootMIP3':    Name,Units = 'Bld 3 IP Root Moment','kN-m'
    elif Parm == 'RootMOoP1':   Name,Units = 'Bld 1 OoP Root Moment','kN-m'
    elif Parm == 'RootMOoP2':   Name,Units = 'Bld 2 OoP Root Moment','kN-m'
    elif Parm == 'RootMOoP3':   Name,Units = 'Bld 3 OoP Root Moment','kN-m'
    elif Parm == 'BldPitch1':   Name,Units = 'Bld 1 Pitch Angle','deg'
    elif Parm == 'BldPitch2':   Name,Units = 'Bld 2 Pitch Angle','deg'
    elif Parm == 'BldPitch3':   Name,Units = 'Bld 3 Pitch Angle','deg'
    elif Parm == 'HSShftTq':    Name,Units = 'High-Speed Shaft Torque','kN-m'
        
    return Name, Units

# =============================================================================
# --------------------------- MISCELLANEOUS -----------------------------------
# =============================================================================

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
    


def Turb2Grid(TurbName):
    """ Load grid parameters for given turbine
    """	
    if TurbName == 'WP0.75A08V00':
        ZHub = 60.
        DZ   = 56.
        DY   = 56.
        n_z  = 15
        n_y  = 15
    elif TurbName == 'WP1.5A08V03':
        ZHub = 84.
        DZ   = 80.
        DY   = 80.
        n_z  = 17
        n_y  = 17

    elif TurbName == 'WP3.0A02V02':
        ZHub = 119.
        DZ   = 120.
        DY   = 120.
        n_z  = 17
        n_y  = 17

    elif TurbName == 'WP5.0A04V00':
        ZHub = 154.
        DZ   = 140.
        DY   = 140.
        n_z  = 15
        n_y  = 15

    else:
        IOError('Turbine {:s} is not coded.'.format(TurbName))

    return (ZHub,DZ,DY,n_z,n_y)          
        
        
def stylepath(style):
    """ Get path to pyplot style file
    """
    
    if (style == 'duke_paper'):
        stylepath = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
                    '2016-02-15_dissertation\\figure_code\\duke_paper.mplstyle'
    elif (style == 'duke_presentation'):
        stylepath = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
                    '2016-02-15_dissertation\\figure_code\\duke_presentation.mplstyle'
    elif (style == 'jmr'):
        stylepath = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
                    '2016-02-15_dissertation\\figure_code\\jmr.mplstyle'
    else:
        ErrStr = 'No defined style type \"{:s}\"'.format(style)
        raise ValueError(ErrStr)

    return stylepath
