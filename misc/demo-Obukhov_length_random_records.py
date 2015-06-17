"""
Demonstration of calculation of ``local'' Obukhov length for random records
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import scipy.io as scio
import random

n_recs  = 1                                 # no. random records to draw
yearRng = [2012,2015]                       # range of years


i_rec = 0
while i_rec < n_recs:
    
    year  = random.randint(yearRng[0],yearRng[1])
    month = random.randint(1,12)
    day   = random.randint(1,31)
    hour  = random.randint(0,23)
    minut = random.randint(0,6)*10
    time_tup = (year,month,day,hour,minut)
    print(time_tup)
    fpath = jr.NRELtime2fpath(time_tup)
    print(fpath)
    
    i_rec += 1
