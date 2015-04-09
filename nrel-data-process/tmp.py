
def dateToM4fname(datetimeS):
    """ Convert array of date information to wildacrd name for M4
        dataset.

        Args:
            datetime (list): [year,month,day,hour,minute], strings

        Returns:
            fname (string): beginning of filename
    """
    
    fname = datetimeS[1] + '_' + datetimeS[2] + '_' + datetimeS[0] \
            + '_' + datetimeS[3] + '_' + datetimeS[4]
    
    return fname


    


# ==============================================================================
import os
import glob
import fnmatch
import scipy.io as sio

# directory path
basedir = '/media/jrinker/JRinker SeaGate External/data/nrel-20Hz/'

# date information
datetime = ['2013', '02', '03', '00', '10']
height = '15'

# get list of files
fnameWC = dateToM4fname(datetime)
dname = os.path.join(basedir, datetime[0], datetime[1], datetime[2])
fname = []
for file in os.listdir(dname):
    if (file.endswith('.mat') and fnameWC in file):
        fname = file
if not fname:
    print 'ERROR: file does not exist'

# load mat file
struc = sio.loadmat(os.path.join(dname,fname))

# extract parameters



