""" Create Python dictionary with metadata information.

    Jenni Rinker, 09-Apr-2015.
"""
def main():
    import os
    import numpy as np
    import time, calendar
    import nrel_lib.nrel_io as nio
    import scipy.io as scio

    # define base directories
    basedir20Hz = '/media/jrinker/JRinker SeaGate External/data/nrel-20Hz/'
    basedir10m  = '/media/jrinker/JRinker SeaGate External/data/nrel-10min/'

    # define useful things
##    fields = nio.metadataFields()
    heights    = np.array([15,30,50,76,100,131])

    # initialize data table
    data = np.empty([0,0])

    # data columns
##    iRecTime = fields.index('Record_Time')

    # recursively loop through all .mat files in 20 Hz directory
    for root, dirs, files in os.walk(basedir20Hz,topdown=False):
        for fname in files:
            if fname.endswith('.mat'):

                # load in 20 Hz datafile
                struc20 = scio.loadmat(os.path.join(root,fname))
                
                # calculate time stamps
                timestamp = nio.fname2time(fname)
                timenow   = calendar.timegm(time.gmtime())

                # loop through heights
                for height in heights:

                    # cup speed
                    nio.interpCupSpeed(struc20,height)

                    # wind direction

                    # precipitation

                    # sonic readings

                    # calculations

                    # create new row

                    # append row
                
                return



if (__name__ == '__main__'):
    main()
