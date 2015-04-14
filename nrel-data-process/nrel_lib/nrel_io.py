""" Functions for NREL data analysis
"""
# ==============================================================================
def fname2time(fname):
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

    # arrange timestamp info in tuple
    time_tup = (year,month,day,hour,minute,0,0,0)

    # convert tuple to float
    time_flt = calendar.timegm(time_tup)

    return time_flt


# ==============================================================================
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
              'up_wp','vp_wp','wp_Tp','up_vp','tau_u','tau_v','tau_w']
    else:
        print '***ERROR***: that dataset is not coded yet.'
        
    return fields


# ==============================================================================
def getBasedir(dataset):
    """ Get path to base directory and check if it exists
    """
    import os

    if (dataset == 'NREL'):
        basedir = '/media/jrinker/JRinker SeaGate External/data/nrel-20Hz/'
        if not os.path.exists(basedir):
            print '***ERROR***: base directory does not exist.'
            del basedir

    else:
        print '***ERROR***: that dataset is not coded yet.'

    return basedir


# ==============================================================================
def makemetadata(dataset):
    """ Construct metadata table
    """

    # get base directory, check it exists
    basedir = getBasedir(dataset)

    # process NREL dataset
    if (dataset == 'NREL'):
        metadata = makeNRELmetadata(basedir)

    else:
        print '***ERROR***: that dataset is not coded yet.'

    return metadata


# ==============================================================================
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
                    row[0,0] = fname2time(fname)
                    row[0,1] = calendar.timegm(time.gmtime())
                    row[0,2] = height
                    metadata = np.vstack([metadata,row])
                
    return metadata


# ==============================================================================
def extractNRELparameters(struc,height):
    """
    """
    import numpy as np

    dataset = 'NREL'

    fields = metadataFields(dataset)

    parameters = np.empty([1,len(fields)])

    for i in range(len(fields)):
        field = fields[i]
        parameters[0,i] = calculatefield(dataset,struc,height,field)
            
    return parameters


# ==============================================================================
def calculatefield(dataset,struc,ht,field):
    """
    """
    import sys
    sys.path.append('/home/jrinker/git/dissertation/')
    import JR_Library.ExtractWindParameters as jr
    import numpy as np

    if (dataset == 'NREL'):
        dt = 0.05

        try:
        
            if ((field == 'Record_Time') or (field == 'Processed_Time') \
                or (field == 'Height')):
                value = 0
                
            elif (field == 'Wind_Speed_Cup'):
                hts_cup = np.array([10,26,80,88,134])
                upperHt = hts_cup[np.where(ht < hts_cup)[0][0]]
                lowerHt = hts_cup[np.where(ht < hts_cup)[0][0] - 1]
                lowerField = 'Cup_WS_' + str(lowerHt) + 'm'
                upperField = 'Cup_WS_' + str(upperHt) + 'm'
                upperWS = np.nanmean(struc[upperField][0,0][0])
                lowerWS = np.nanmean(struc[lowerField][0,0][0])
                value = np.interp(ht,[lowerHt,upperHt], \
                                  [lowerWS,upperWS])
                
            elif (field == 'Wind_Direction'):
                hts_vane = np.array([10,26,88,134])
                upperHt = hts_vane[np.where(ht < hts_vane)[0][0]]
                lowerHt = hts_vane[np.where(ht < hts_vane)[0][0] - 1]
                lowerField = 'Vane_WD_' + str(lowerHt) + 'm'
                upperField = 'Vane_WD_' + str(upperHt) + 'm'
                lowerWD = struc[upperField][0,0][0]
                upperWD = struc[lowerField][0,0][0]
                lowerWD_mean = np.angle(np.nanmean(np.exp(1j*lowerWD \
                            *np.pi/180)))*180/np.pi
                upperWD_mean = np.angle(np.nanmean(np.exp(1j*upperWD \
                            *np.pi/180)))*180/np.pi
                value = np.mod(np.angle(np.nanmean(np.exp(1j* \
                    np.array([lowerWD_mean, upperWD_mean])* \
                    np.pi/180)))*180/np.pi,360)

            elif (field == 'Precipitation'):
                value = np.nanmean(struc['PRECIP_INTEN'][0,0][0])

            elif (field == 'Mean_Wind_Speed'):
                value = np.nanmean(struc['Sonic_u_' + str(ht) + 'm'][0,0][0])

            elif (field == 'Sigma_u'):
                value = np.nanstd(struc['Sonic_u_' + str(ht) + 'm'][0,0][0])

            elif (field == 'Concentration_u'):
                u = struc['Sonic_u_' + str(ht) + 'm'][0,0][0]
                rho, mu = jr.signalPhaseCoherence(u)
                value = rho

            elif (field == 'Location_u'):
                u = struc['Sonic_u_' + str(ht) + 'm'][0,0][0]
                rho, mu = jr.signalPhaseCoherence(u)
                value = mu

            elif (field == 'Sigma_v'):
                value = np.nanstd(struc['Sonic_v_' + str(ht) + 'm'][0,0][0])

            elif (field == 'Concentration_v'):
                v = struc['Sonic_v_' + str(ht) + 'm'][0,0][0]
                rho, mu = jr.signalPhaseCoherence(v)
                value = rho

            elif (field == 'Location_v'):
                v = struc['Sonic_v_' + str(ht) + 'm'][0,0][0]
                rho, mu = jr.signalPhaseCoherence(v)
                value = mu

            elif (field == 'Sigma_w'):
                value = np.nanstd(struc['Sonic_w_' + str(ht) + 'm'][0,0][0])

            elif (field == 'Concentration_w'):
                w = struc['Sonic_w_' + str(ht) + 'm'][0,0][0]
                rho, mu = jr.signalPhaseCoherence(w)
                value = rho

            elif (field == 'Location_w'):
                w = struc['Sonic_w_' + str(ht) + 'm'][0,0][0]
                rho, mu = jr.signalPhaseCoherence(w)
                value = mu

            elif (field == 'up_wp'):
                u = struc['Sonic_u_' + str(ht) + 'm'][0,0][0]
                w = struc['Sonic_w_' + str(ht) + 'm'][0,0][0]
                value = np.nanmean((u-np.nanmean(u))*(w-np.nanmean(w)))           

            elif (field == 'vp_wp'):
                v = struc['Sonic_v_' + str(ht) + 'm'][0,0][0]
                w = struc['Sonic_w_' + str(ht) + 'm'][0,0][0]
                value = np.nanmean((v-np.nanmean(v))*(w-np.nanmean(w)))

            elif (field == 'wp_Tp'):
                w = struc['Sonic_w_' + str(ht) + 'm'][0,0][0]
                T = struc['Sonic_Temp_rotated_' + str(ht) + 'm'][0,0][0]
                value = np.nanmean((w-np.nanmean(w))*(T-np.nanmean(T)))

            elif (field == 'up_vp'):
                u = struc['Sonic_u_' + str(ht) + 'm'][0,0][0]
                v = struc['Sonic_v_' + str(ht) + 'm'][0,0][0]
                value = np.nanmean((u-np.nanmean(u))*(v-np.nanmean(v))) 

            elif (field == 'tau_u'):
                u = struc['Sonic_u_' + str(ht) + 'm'][0,0][0]
                value = jr.calculateKaimal(u,dt)

            elif (field == 'tau_v'):
                v = struc['Sonic_v_' + str(ht) + 'm'][0,0][0]
                value = jr.calculateKaimal(v,dt)

            elif (field == 'tau_w'):
                w = struc['Sonic_w_' + str(ht) + 'm'][0,0][0]
                value = jr.calculateKaimal(w,dt)

            else:
                print '***ERROR***: field {} is not coded for \"NREL\".'.format(field)

        except KeyError:
            print '***WARNING***: KeyError for {}'.format(field)
            value = float('nan')
        
    else:
        print '***ERROR***: that dataset is not coded yet.'

    return value


# ==============================================================================
def loadMetadata(fname):
    """
    """
    import numpy as np

    metadata = np.loadtxt(fname,delimiter=',',skiprows=1)

    with open(fname,'r') as f:
        header = f.readline()

    fields = header.lstrip('# ').rstrip('\n').split(',')

    return (fields,metadata)

    

