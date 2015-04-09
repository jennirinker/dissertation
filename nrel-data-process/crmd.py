""" Create Python dictionary with metadata information.

    Jenni Rinker, 09-Apr-2015.
"""
import os
import numpy as np
from nrel_lib.nrel_io import *

# define base directory
basedir = '/media/jrinker/JRinker SeaGate External/data/nrel-20Hz/'

# define metadata fields
fields = metadataFields()

# initialize data table
data = np.empty([0,0])

# data columns
iRecTime = fields.index('Record_Time')

# for each year
year = 2015

    # for each month
month = 02

        # for each day
day = 06

            # get list of records for that day
dname = makeDataPath(basedir,year,month,day)
matFiles = listMatFiles(dname)

            # for each record
irec = 0

fname = matFiles[irec]

                # calculate the floating time-stamp
timestamp = fname2time(fname)

                # if the record does not exist in the datatable


                    # extract the parameters
