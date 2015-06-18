"""
Demonstration of calculation of ``local'' Obukhov length for random records
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import scipy.io as scio
import random
import calendar, time
import numpy as np


n_recs  = 1                                 # no. random records to draw
yearRng = [2012,2015]                       # range of years
fields = jr.metadataFields('NREL')          # MD fields

CSLim  = 3                          # lower cup speed limit
dir1   = 240                        # CCW edge for direction range
dir2   = 315                        # CW edge for direction range
preLim = 2.7                        # lower precipitation limit

# column indices for each value
CScol  = fields.index('Wind_Speed_Cup')
dirCol = fields.index('Wind_Direction')
preCol = fields.index('Precipitation')

i_rec = 0
metadata = np.empty((n_recs,len(fields)))
while i_rec < n_recs:
    
    # draw random time stamp/height
    year   = random.randint(yearRng[0],yearRng[1])
    month  = random.randint(1,12)
    day    = random.randint(1,31)
    hour   = random.randint(0,23)
    minute = random.randint(0,6)*10
    height = random.choice([15,30,50,76,100,131])
    time_tup = (year,month,day,hour,minute,height)
        
    # check if it exists
    fpath = jr.NRELtime2fpath(time_tup)
    if len(fpath) > 0:
        
        # load structure
        struc = scio.loadmat(fpath)
        
        # calculate parameters
        row = jr.extractNRELparameters(struc,height)
        row[0] = jr.timetup2flt(time_tup)
        row[1] = calendar.timegm(time.gmtime())
        row[2] = height
        
        flagCS   = row[CScol] > CSLim
        flagDir1 = row[dirCol] >= dir1
        flagDir2 = row[dirCol] <= dir2
        flagPrec = row[preCol] >= preLim
        
        if (flagCS and flagDir1 and flagDir2 and flagPrec):
            metadata[i_rec] = row.reshape(1,len(fields))
    
            i_rec += 1
