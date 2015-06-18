"""
Library of functions used in dissertation analyses.

Module is separated into groups of functions:
    - File I/O:
        Loading metadata files, TurbSim files, data records
    - Metadata processing:
        Extracting wind parameters from data
    - Data analysis:
        Calculating NSAE,
        fitting composite distributions
    - Simulation:
        Kaimal simulation
    - IEC:
        Kaimal spectrum, spatial coherence, etc.
    - Miscellaneous:
        Other functions

Jenni Rinker, Duke University
"""

# ==============================================================================
# FILE I/O
# ==============================================================================

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

    # if it is a text file
    if (fname.endswith('txt')):
    
        # strip fields from header line
        with open(fname,'r') as f:
            fields = f.readline().lstrip('# ').rstrip('\n').split(',')
            
        # read in data table
        metadata = np.loadtxt(fname,delimiter=',',skiprows=1)

    # if it is a .mat file
    elif (fname.endswith('mat')):

        # code this later
        print('***ERROR*** Have not yet coded loading a .mat file')
    
    return (fields, metadata)


def loadNRELmatlab():
    """ Load the metadata table processed in Matlab

        Returns:
            fields (list): names of fields in coumns
            metadata (numpy array): array of atmospheric params
    """
    import scipy.io as scio
    import numpy as np
    
    # path to matlab-processed metadata table
    matpath = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                  'processed_data\\NREL_metadata_mat.mat'

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
    fields[fields.index('Cup_speed')] = 'Wind_Speed_Cup'
    fields[fields.index('Wind_direction')] = 'Wind_Direction'
    fields[fields.index('Precip_intens')] = 'Precipitation'
    fields[fields.index('rho_u')] = 'Concentration_u'

    return (fields, metadata)

def loadtimeseries(dataset,timestamp,ht):
    """ Load the turbulent time series for a given dataset and timestamp.

        Args:
            dataset (string): flag to indicate dataset
            timestamp (float/tuple): float or tuple representing time value

        Returns:
            t (numpy array): array of time values
            u (numpy array): array of longitudinal win velocities
    """
    import scipy.io as scio
    import numpy as np
    import calendar

    if (dataset == 'NREL'):

        # time parameters
        N  = 12000
        dt = 0.05

        # convert tuple to timestamp if necessary
        if (type(timestamp) == tuple):
            timestamp = timetup2flt(timestamp)

        # skip to last step if it's a structure
        if (type(timestamp) == dict):
            struc = timestamp
        else:
            # get file path
            fpath = NRELtime2fpath(timestamp)
            struc = scio.loadmat(fpath)

        # load data
        t = np.arange(12000)*dt
        u = np.squeeze(struc['Sonic_u_' + str(ht) + 'm'][0,0][0])
        v = np.squeeze(struc['Sonic_v_' + str(ht) + 'm'][0,0][0])
        w = np.squeeze(struc['Sonic_w_' + str(ht) + 'm'][0,0][0])


    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return (t,u, v, w)


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

# ==============================================================================
# METADATA PROCESSING
# ==============================================================================

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
        
        # filter out the rows with NaN values???        
        
        # screen remaining data
        cleandata = metadata[np.where(metadata[:,CScol] > CSLim)]
        cleandata = cleandata[np.where(cleandata[:,dirCol] >= dir1)]
        cleandata = cleandata[np.where(cleandata[:,dirCol] <= dir2)]
        cleandata = cleandata[np.where(cleandata[:,preCol] >= preLim)]
        
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)
        
    return cleandata

def NRELfname2time(fname):
    """ Convert filename to float representing time stamp

        Args:
            fname (string): name of file

        Returns:
            timestamp (float): float representing timestamp
    """
    import calendar

    # extract date information
    year   = int(fname[6:10])
    month  = int(fname[0:2])
    day    = int(fname[3:5])
    hour   = int(fname[11:13])
    minute = int(fname[14:16])

    # convert tuple to float
    time_flt = timetup2flt(time_tup)

    return time_flt


def NRELtime2fpath(time_flt):
    """ Convert float to path to data structure

        Args:
            time_flt (float): float representing timestamp

        Returns:
            fpath (string): path to data 
            
    """
    import time
    import os
    import glob

    # check if it's a tuple
    if type(time_flt) == tuple:
        time_tup = time_flt
    else:
        time_tup = timeflt2tup(time_flt)

    # convert tuple to string values
    yearS  = str(time_tup[0])
    monthS = str(time_tup[1]).zfill(2)
    dayS   = str(time_tup[2]).zfill(2)
    hourS  = str(time_tup[3]).zfill(2)
    minS   = str(time_tup[4]).zfill(2)

    # get directory path
    basedir = getBasedir('NREL')
    dirpath = os.path.join(basedir,yearS,monthS,dayS)    
    # get filename
    fname_part = '_'.join([monthS,dayS,yearS,hourS,minS])
    fpath = glob.glob(os.path.join(dirpath,fname_part)+'*.mat')

    # check if file doesn't exist
    if len(fpath) > 0:  fpath = fpath[0]
    else:               fpath = ''
    
    return fpath


def metadataFields(dataset):
    """ Define list of fields to be stored in metadata table

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
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)
        
    return fields


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


def makeNRELmetadata(basedir):
    """ Make metadata table for NREL data
    """
    import numpy as np
    import scipy.io as scio
    import os
    import calendar, time

    # array of sampling heights
    heights = np.array([15,30,50,76,100,131])

    # initialize array
    metadata = np.empty([0,len(metadataFields('NREL'))])

    # recursively loop through all .mat files in 20 Hz directory
    for root, dirs, files in os.walk(basedir,topdown=False):
        print ' Procesing ' + root
        for fname in files:
            if fname.endswith('.mat'):

                # load in 20 Hz datafile
                struc = scio.loadmat(os.path.join(root,fname))                
                # loop through heights
                for height in heights:

                    row = extractNRELparameters(struc,height)
                    row[0] = NRELfname2time(fname)
                    row[1] = calendar.timegm(time.gmtime())
                    row[2] = height
                    metadata = np.vstack([metadata,row])
                
    return metadata


def extractNRELparameters(struc,height):
    """
    """
    import numpy as np

    dataset = 'NREL'

    fields = metadataFields(dataset)

    parameters = np.empty(len(fields))

    outdict = calculatefield(dataset,struc,height)

    for i in range(len(fields)):
        field = fields[i]
        parameters[i] = outdict[field]
            
    return parameters


def interpolateparameter(dataset,struc,ht,field):
    """ Interpolate a parameter value

        Args:
            dataset (string): flag to indicate which dataset to analyze
            struc (numpy structure): loaded from matlab
            ht (float): measurement height to interpolate
            field (string): field of interest for interpolation

        Returns:
            value (float): interpolated value
    """
    import numpy as np

    if (dataset == 'NREL'):
        dt = 0.05

        try:

            if (field == 'Wind_Speed_Cup'):

                # define measurement heights
                msmnt_hts = np.array([10,26,80,88,134])
                fld_str   = 'Cup_WS_'

                # force height to be in range
                ht = min(msmnt_hts[-1],ht)
                ht = max(msmnt_hts[0],ht)

                # perform interpolation
                upperHt = msmnt_hts[np.where(ht < msmnt_hts)[0][0]]
                lowerHt = msmnt_hts[np.where(ht < msmnt_hts)[0][0] - 1]
                lowerField = fld_str + str(lowerHt) + 'm'
                upperField = fld_str + str(upperHt) + 'm'
                upperVal   = np.nanmean(struc[upperField][0,0][0])
                lowerVal   = np.nanmean(struc[lowerField][0,0][0])
                value = np.interp(ht,[lowerHt,upperHt], \
                                    [lowerVal,upperVal])

            elif (field == 'Temperature'):
                
                # define measurement heights
                msmnt_hts = np.array([3,26,88])
                fld_str   = 'Air_Temp_'

                # force height to be in range
                ht = min(msmnt_hts[-1],ht)
                ht = max(msmnt_hts[0],ht)

                # perform interpolation
                upperHt = msmnt_hts[np.where(ht <= msmnt_hts)[0][0]]
                lowerHt = msmnt_hts[np.where(ht <= msmnt_hts)[0][0] - 1]
                lowerField = fld_str + str(lowerHt) + 'm'
                upperField = fld_str + str(upperHt) + 'm'
                upperVal   = np.nanmean(struc[upperField][0,0][0])
                lowerVal   = np.nanmean(struc[lowerField][0,0][0])
                value = np.interp(ht,[lowerHt,upperHt], \
                                    [lowerVal,upperVal])
                
            elif (field == 'Wind_Direction'):
                
                # define measurement heights
                msmnt_hts = np.array([10,26,88,134])
                fld_str   = 'Vane_WD_'

                # force height to be in range
                ht = min(msmnt_hts[-1],ht)
                ht = max(msmnt_hts[0],ht)

                # perform interpolation
                upperHt = msmnt_hts[np.where(ht < msmnt_hts)[0][0]]
                lowerHt = msmnt_hts[np.where(ht < msmnt_hts)[0][0] - 1]
                lowerField = fld_str + str(lowerHt) + 'm'
                upperField = fld_str + str(upperHt) + 'm'
                lowerVal = struc[upperField][0,0][0]
                upperVal = struc[lowerField][0,0][0]
                lower_mean = np.angle(np.nanmean(np.exp(1j*lowerVal \
                            *np.pi/180)))*180/np.pi
                upper_mean = np.angle(np.nanmean(np.exp(1j*upperVal \
                            *np.pi/180)))*180/np.pi
                value = np.mod(np.angle(np.nanmean(np.exp(1j* \
                    np.array([lower_mean, upper_mean])* \
                    np.pi/180)))*180/np.pi,360)

            else:
                errStr = 'Field \"{}\" is not coded for dataset \"NREL\".'.format(field)
                raise AttributeError(errStr)            

        except KeyError:
            print '***WARNING***: KeyError for {}'.format(field)
            value = float('nan')
        
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return value
        
def closestparameter(dataset,struc,ht,field):
    """ Return the closest parameter value

        Args:
            dataset (string): flag to indicate which dataset to analyze
            struc (numpy structure): loaded from matlab
            ht (float): measurement height to interpolate
            field (string): field of interest for interpolation

        Returns:
            value (float): closest value
    """
    import numpy as np

    if (dataset == 'NREL'):
        dt = 0.05

        try:

            if (field == 'Temperature'):
                
                # define measurement heights
                msmnt_hts = np.array([3,26,88])
                fld_str   = 'Air_Temp_'

                # force height to be in range
                ht = min(msmnt_hts[-1],ht)
                ht = max(msmnt_hts[0],ht)

                # find closest value
                upperHt = msmnt_hts[np.where(ht <= msmnt_hts)[0][0]]
                lowerHt = msmnt_hts[np.where(ht <= msmnt_hts)[0][0] - 1]
                dh      = upperHt - lowerHt
                ht_clst = int(round((ht-lowerHt)/dh)*dh + lowerHt)
                clstField = fld_str + str(lowerHt) + 'm'
                value   = np.nanmean(struc[clstField][0,0][0])

            else:
                errStr = 'Field \"{}\" is not coded for dataset \"NREL\".'.format(field)
                raise AttributeError(errStr)            

        except KeyError:
            print '***WARNING***: KeyError for {}'.format(field)
            value = float('nan')
        
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return value
        

def calculatefield(dataset,struc,ht):
    """ Save atmophseric parameters in output dictionary

    """
    import sys
    import numpy as np

    if (dataset == 'NREL'):
        dt = 0.05

        # calculate fields necessary for later calculations
        (t,u,v,w) = loadtimeseries(dataset,struc,ht)
        T_sonic_K = np.squeeze(struc['Sonic_Temp_rotated_' + \
                             str(ht) + 'm'][0,0][0]) + 298.15
        up        = nandetrend(t,u)
        vp        = nandetrend(t,v)
        wp        = nandetrend(t,w)
        Tp        = nandetrend(t,T_sonic_K)
        upwp_bar  = np.nanmean(up*wp)
        upvp_bar  = np.nanmean(up*vp)
        vpwp_bar  = np.nanmean(vp*wp)
        wpTp_bar  = np.nanmean(wp*Tp)
        ustar     = ( (upwp_bar)**2 + (vpwp_bar)**2 ) ** 0.25
        rhou,muu  = signalPhaseCoherence(up)
        rhov,muv  = signalPhaseCoherence(vp)
        rhow,muw  = signalPhaseCoherence(wp)
        Tbar_in_K = interpolateparameter(dataset, \
                            struc,ht,'Temperature') + 298.15
        Tbar_ne_K = closestparameter(dataset, \
                            struc,ht,'Temperature') + 298.15

        # initialize output dictionary
        outdict = {}

        # save values
        outdict['Record_Time']     = 0
        outdict['Processed_Time']  = 0
        outdict['Height']          = 0
        outdict['Wind_Speed_Cup']  = interpolateparameter( \
                                        dataset,struc,ht,'Wind_Speed_Cup')
        outdict['Wind_Direction']  = interpolateparameter( \
                                        dataset,struc,ht,'Wind_Direction')
        outdict['Precipitation']   = np.nanmean( \
                                        struc['PRECIP_INTEN'][0,0][0])
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
                                            np.mean(u),dt)
        outdict['tau_v']           = calculateKaimal(vp,dt)
        outdict['tau_w']           = calculateKaimal(wp,dt)
        outdict['MO_Length_interp'] = -(Tbar_in_K * ustar**3) \
                                     /(0.41 * 9.81 * wpTp_bar)
        outdict['MO_Length_near']  = - (Tbar_ne_K * ustar**3)/ \
                                     (0.41 * 9.81 * wpTp_bar)
        
    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return outdict

##def calculatefield(dataset,struc,ht,field):
##    """
##    """
##    import sys
##    import numpy as np
##
##    if (dataset == 'NREL'):
##        dt = 0.05
##
##        (t,u,v,w) = loadtimeseries(dataset,struc,ht)
##        up        = nandetrend(t,u)
##        vp        = nandetrend(t,v)
##        wp        = nandetrend(t,w)
##
##        try:
##        
##            if ((field == 'Record_Time') or (field == 'Processed_Time') \
##                or (field == 'Height')):
##                value = 0
##                
##            elif (field == 'Wind_Speed_Cup'):
##                value = interpolateparameter(dataset,struc,ht,field)
##                
##            elif (field == 'Wind_Direction'):
##                value = interpolateparameter(dataset,struc,ht,field)
##
##            elif (field == 'Precipitation'):
##                value = np.nanmean(struc['PRECIP_INTEN'][0,0][0])
##
##            elif (field == 'Mean_Wind_Speed'):
##                value = np.nanmean(u)
##
##            elif (field == 'Sigma_u'):
##                value = np.nanstd(up)
##
##            elif (field == 'Concentration_u'):
##                value = signalPhaseCoherence(up)[0]
##
##            elif (field == 'Location_u'):
##                value = signalPhaseCoherence(up)[1]
##
##            elif (field == 'Sigma_v'):
##                value = np.nanstd(vp)
##
##            elif (field == 'Concentration_v'):
##                value = signalPhaseCoherence(vp)[0]
##
##            elif (field == 'Location_v'):
##                value = signalPhaseCoherence(vp)[1]
##
##            elif (field == 'Sigma_w'):
##                value = np.nanstd(wp)
##
##            elif (field == 'Concentration_w'):
##                value = signalPhaseCoherence(wp)[0]
##                
##            elif (field == 'Location_w'):
##                value = signalPhaseCoherence(wp)[1]
##
##            elif (field == 'up_wp'):
##                value = np.nanmean(up*wp)           
##
##            elif (field == 'vp_wp'):
##                value = np.nanmean(vp*wp)
##
##            elif (field == 'wp_Tp'):
##                T = np.squeeze(struc['Sonic_Temp_rotated_' + \
##                                     str(ht) + 'm'][0,0][0]) + 298.15
##                Tp = nandetrend(t,T)
##                value = np.nanmean(wp*Tp)
##
##            elif (field == 'up_vp'):
##                value = np.nanmean(up*vp) 
##
##            elif (field == 'tau_u'):
##                value = calculateKaimal(up + np.mean(u),dt)
##
##            elif (field == 'tau_v'):
##                value = calculateKaimal(vp,dt)
##
##            elif (field == 'tau_w'):
##                value = calculateKaimal(wp,dt)
##
##            elif (field == 'MO_Length_interp'):
##                T_bar = interpolateparameter(dataset, \
##                    struc,ht,'Temperature') + 298.15
##                T     = np.squeeze(struc['Sonic_Temp_rotated_' + \
##                          str(ht) + 'm'][0,0][0]) + 298.15
##                Tp       = nandetrend(t,T)
##                upwp_bar = np.nanmean(up*wp)
##                vpwp_bar = np.nanmean(vp*wp)
##                wpTp_bar = np.nanmean(wp*Tp)
##                ustar    = ( (upwp_bar)**2 + (vpwp_bar)**2 ) ** 0.25
##                value    = - (T_bar * ustar**3)/(0.41 * 9.81 * wpTp_bar)
##
##            elif (field == 'MO_Length_closest'):
##                T_bar = closestparameter(dataset, \
##                    struc,ht,'Temperature') + 298.15
##                T     = np.squeeze(struc['Sonic_Temp_rotated_' + \
##                            str(ht) + 'm'][0,0][0]) + 298.15
##                Tp       = nandetrend(t,T)
##                upwp_bar = np.nanmean(up*wp)
##                vpwp_bar = np.nanmean(vp*wp)
##                wpTp_bar = np.nanmean(wp*Tp)
##                ustar    = ( (upwp_bar)**2 + (vpwp_bar)**2 ) ** 0.25
##                value    = - (T_bar * ustar**3)/(0.41 * 9.81 * wpTp_bar)
##                
##            else:
##                errSt = 'Field {} is not coded for ' + \
##                        'dataset \"NREL\".'.format(field)
##                raise AttributeError(errStr)
##
##        except KeyError:
##            print '***WARNING***: KeyError for {}'.format(field)
##            value = float('nan')
##        
##    else:
##        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
##        raise AttributeError(errStr)
##
##    return value


# ==============================================================================
# DATA ANALYSIS
# ==============================================================================

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

def samplePhaseCoherence(theta):
    """ Return concentration and location parameters for a sample of wrapping
        random variables.
        
        Args:
            theta (numpy array): 1D array of sample of angles
       
        Returns:
            rho (float): concentration parameter
            mu (float): location parameter
    """
    import numpy as np
    
    # if (2+)D array is fed in, halt with error
    if ((len(theta.shape)>1) and (theta.shape[0] != 1 and theta.shape[1] != 1)):
        print 'ERROR: samplePhaseCoherence only works on 1D arrays'
        return []
    
    n_t = theta.size                        # number of elements in sample
    V = np.sum(np.exp(1j*theta))/n_t        # mean resultant vector
    rho = np.abs(V)                         # concentration parameter
    mu  = np.angle(V)                       # location parameter
    
    return (rho,mu) 


def signalPhaseDifferences(x):
    """ Return phase differences for a time history.
    
        Args:
            x (numpy array): 1D array of time history
       
        Returns:
            rho (float): concentration parameter
            mu (float): location parameter
    """
    import numpy as np
    
    # if (2+)D array is fed in, halt with error
    if ((len(x.shape)>1) and (x.shape[0] != 1 and x.shape[1] != 1)):
        errStr = 'signalPhaseDifferences only works on 1D arrays'
        raise AttributeError(errStr)
    
    n_t = x.shape[0]                        # no. of total components
    n_f = uniqueComponents(n_t)             # no. of unique components

    X = np.fft.fft(x)/n_t                   # Fourier vector
    Xuniq = X[:n_f]                         # unique Fourier components
    dtheta = np.angle(np.divide( \
        Xuniq[1:],Xuniq[:-1]))              # phase differences
    
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
        rho, mu = samplePhaseCoherence( \
            dtheta)                                 # temp coh parameters
    
    return (rho,mu)


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


def spatialCoherence(Xi,Xj,df):
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
        print 'ERROR: calculateKaimal only works on 1D arrays'
        return []

    # grid search parameters
    n_m = 200;                                  # number of points in grid
    tau_l = -1;                                 # left tau coefficient
    tau_r = 3;                                  # right tau coefficient

    # intermediate parameters
    n_t  = x.shape[0]                           # no. of time steps
    T    = n_t*dt                               # total time
    n_f  = uniqueComponents(n_t)                # no. unique components
    df   = 1/T                                  # frequency resolution
    f    = df*np.arange(0,n_f). \
        reshape(n_f,1)                          # frequency vector
    sig  = np.std(x,ddof=1)                     # std. deviation

    X_dat  = np.fft.fft(x,axis=0). \
        reshape(n_t,1)/n_t                      # Fourier vector of data
    Sk_dat = X2Sk(X_dat)[1:]* \
        np.ones((1,n_m))                        # data discrete PSD from f1 up

    taugrid = np.logspace(tau_l,tau_r,n_m). \
        reshape(1,n_m)                          # grid of tau for search
    
    S_kaim = KaimalSpectrum(f,taugrid,sig)      # continuous Kaimal spectrum
    Sk_kaim = S_kaim * df                       # discrete Kaimal spectrum
    Sk_theo = np.empty(Sk_kaim.shape)
    for i in range(n_m):
        alpha = spectralScale(Sk_kaim[:,i],sig,n_t)
        Sk_theo[:,i] = (alpha**2)*Sk_kaim[:,i]
    Sk_theo = Sk_theo[1:]                       # strip off DC component
    
    J = np.sum(np.power(Sk_dat-Sk_theo, \
        2),axis=0)                              # array of squared errors

    i_min = np.argmin(J)                        # index of min error

    tau = taugrid[0,i_min]                      # tau_minError
    
    return tau


def calculateTurbSimSC(fname,rsep):
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


def getBasedir(dataset):
    """ Get path to base directory and check if it exists
    """
    import os
    import platform

    if (dataset == 'NREL'):
        if (platform.system() == 'Linux'):
            basedir = '/media/jrinker/JRinker SeaGate External/data/nrel-20Hz/'
        elif (platform.system() == 'Windows'):
            basedir = 'G:\\data\\nrel-20Hz'
        if not os.path.exists(basedir):
            errStr = 'Incorrect or unavailable base ' + \
                     'directory for dataset \"{}\".'.format(dataset)
            raise IOError(errStr)

    else:
        errStr = 'Dataset \"{}\" is not coded yet.'.format(dataset)
        raise AttributeError(errStr)

    return basedir


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
        x = np.ones(n_t,n_m)*U;
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

# ==============================================================================
# IEC
# ==============================================================================

def VelProfile(z,zhub,Vhub):
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


def IrefFromClass(turbc):
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


def calcSigma1(Iref,Vhub):
    """ Longitudinal standard deviation
    
        Args:
            Iref (float): reference hub-height turbulence intensity
            
        Returns:
            sigma1 (float): longitudinal standard deviation
    """
    
    sigma1 = Iref*(0.75 * Vhub + 5.6);
    
    return sigma1


def TiProfile(z,zhub,Vhub,turbc):
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


def calcLambda1(zhub):
    """ Longitudinal scale parameter
    
        Args:
            zhub (float): hub-height of turbine in meters
        
        Returns:
            Lambda1 (float): longitudinal scale parameter
    """
    
    Lambda1 = (0.7 * zhub)*(zhub <= 60) \
        + 42*(zhub > 60);               # longitudinal scale parameter
    
    return Lambda1


def SpatialCoherence(zhub,Vhub,rsep,f):
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

    Lambda1 = calcLambda1(zhub)         # longitudinal scale parameter
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


def PSDs(zhub,Vhub,turbc,f):
    """ Continuous power spectral densities for all 3 turbulence components
    
        Args:
            zhub (float): hub-height of turbine in meters
            Vhub (float): hub-height mean velocity in m/s
            turbc (string): turbulence class from IEC 61400-1 Ed. 3
            f (numpy array): frequencies in Hz
            
        Returns:
            Su (numpy array): longitudinal PSD
            Sv (numpy array): lateral PSD
            Sw (numpy array): vertical PSD
    """

    Iref = IrefFromClass(turbc);        # referene turbulence intensity
    sigma1 = calcSigma1(Iref,Vhub);     # longitudinal standard deviation
    sigma2 = 0.8*sigma1;                # lateral standard deviaiton
    sigma3 = 0.5*sigma1;                # vertical standard deviation
    Lambda1 = calcLambda1(zhub);        # longitudinal scale parameter
    L1 = 8.1*Lambda1;                   # longitudinal integral scale   
    L2 = 2.7*Lambda1;                   # lateral integral scale  
    L3 = 0.66*Lambda1;                  # vertical integral scale

    Su = KaimalSpectrum(f,L1/Vhub,sigma1);  # longitudinal spectrum
    Sv = KaimalSpectrum(f,L2/Vhub,sigma2);  # lateral spectrum
    Sw = KaimalSpectrum(f,L3/Vhub,sigma3);  # vertical spectrum

    return (Su,Sv,Sw)

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
    
    theta_out = np.multiply(theta<=np.pi,theta) + \
        np.multiply(theta>np.pi,theta-2*np.pi)  # subtract 2pi from [pi,2pi)
    
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


def timeflt2tup(time_flt):
    """ Convert time in float to tuple

        Args:
            time_flt (float): timestamp in float format

        Returns:
            time_tup (tuple): (year,month,day,hour,minute)
    """
    import time

    time_tup = time.gmtime(time_flt)[:5]

    return time_tup

def timetup2flt(time_tup):
    """ Convert time in float to tuple

        Args:
            time_tup (tuple): (year,month,day,hour,minute)
            
        Returns:
            time_flt (float): timestamp in float format
            
    """
    import calendar

    time_tup += (0,0,0)
    time_flt = calendar.timegm(time_tup)

    return time_flt


def nandetrend(x,y):
    """ Linear detrend ignoring nan values in y

        Args:
            x (numpy array): x values
            y (numpy array): y values

        Returns:
            y_det (numpy array): detrended values with nans replaced
    """
    import numpy as np

    # initialize output
    y_det   = np.empty(y.shape)

    # check if all NaNs
    if np.all(np.isnan(y)):
        y_det[:] = np.NAN

    # linear detrend if 1+ non-nan values
    else:

        # find NaN indices
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

    
##def numpy2latex(a,toprow=None,firstcol=None):
##    """ Print numpy array in LaTeX-friendly formatting
##
##        Args:
##            a (numpy array): array of values to print
##            toprow (list/numpy array): headers for top row
##            firstcol (list/numpy array): headers for first column
##
##        Returns:
##            output (list): list of strings for each row
##    """
##
##    m,n = a.shape                               # get dimensions
##
##    if toprow == None: toprow = ['']*n
##    if firstcol == None: firstcol = ['']*m
##
##    for iRow in range(m):
##
##        rowStr = ''
##
##        for iCol in range(n):

            
                
    

