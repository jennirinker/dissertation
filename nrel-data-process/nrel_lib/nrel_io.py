""" Functions for NREL data analysis
"""

def listMatFiles(dname):
    """ List all .mat files at directory dname

        Args:
            dname (string): path to directory

        Returns:
            matList (list): list of .mat filenames
    """
    import os

    matList = []
    for file in os.listdir(dname):
        if (file.endswith('.mat')):
            matList.append(file)

    return matList

def makeDataPath(basedir,year,month,day):
    """ Construct data path to subfolder

        Ars:
            basedir (string): path to base directory
            year (integer): year of data acquisition
            month (integer): month of data acquisition
            day (integer): day of data acquisition

        Returns:
            dname (string): path to subfolder
    """
    import os

    yearS  = str(year)
    monthS = str(month).zfill(2)
    dayS   = str(day).zfill(2)

    dname = os.path.join(basedir,yearS,monthS,dayS)

    return dname

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

def metadataFields():
    """ Define list of fields to be stored in metadata table

        Returns:
            fields (list): list of strings defining metadata
                columns
    """

    fields = ['Record_Time','Processed_Time','Height','Wind_Speed_Cup', \
              'Wind_Direction','Precipitation','Temperature', \
              'Dewpoint_Temp','Surface_Virtual_Temp','Surface_Pressure', \
              'Mean_Wind_Speed', \
              'Sigma_u','Concentration_u','Location_u', \
              'Sigma_v','Concentration_v','Location_v', \
              'Sigma_w','Concentration_w','Location_w', \
              'up_wp','vp_wp','wp_Tp','up_vp','tau_u','tau_v','tau_w', \
              'Andy_Lu','Andy_Lv','Andy_Lw','Ri_grad_26_88_134m', \
              'Ri_grad_3_10_26_88_134m','Ri_grad_3_10_26_88m', \
              'Ri_grad_10_26_88_134m','Ri_WS_26_88_134m', \
              'Ri_WS_3_10_26_88_134m','Ri_WS_3_10_26_88m', \
              'Ri_WS_10_26_88_134m']

    return fields

def field2varname(field,height):
    """ Convert fieldname to variable name

        Args:
            field (string): field name stored in Python dictionary
            height (integer): measurement height

        Returns:
            varname (string): variable name in data
    """
    htS = str(height)
    
    if (field == 'Wind_Speed_Cup'):
        varname = 'Wind_Speed_Cup_' + htS + 'm'
        
    elif (field == 'Wind_Direction'):
        varname = 'Wind_Direction_Vane_' + htS + 'm_mean'
        
    elif (field == 'Precipitation'):
        varname = 'Raw_PRECIP_INTEN_mean'
        
    elif (field == 'Temperature'):
        varname = 'Air_Temperature_' + htS + 'm'
        
    elif (field == 'Surface_Pressure'):
        varname = 'Raw_Baro_Presr_3m_mean'

    else:
        print 'ERROR: that fieldname is not listed in this function'
        return []

# TODO: add all fieldnames that are used in any calculations (e.g., sonics, etc.)

    return varname
