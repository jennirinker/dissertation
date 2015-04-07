import os
import glob
import fnmatch
import sys
sys.path.append('..')
import JR_Library.data_io as dio
import scipy.io as sio

# directory path
basedir = '/media/jrinker/JRinker SeaGate External/data/nrel-20Hz/'

# date information
datetime = ['2013', '02', '03', '00', '10']

# get list of files
fnameWC = dio.dateToM4fname(datetime)
dname = os.path.join(basedir, datetime[0], datetime[1], datetime[2])
fname = []
for file in os.listdir(dname):
    if (file.endswith('.mat') and fnameWC in file):
        fname = file
if not fname:
    print 'ERROR: file does not exist'

# load mat file
struc = sio.loadmat(os.path.join(dname,fname))
struc
