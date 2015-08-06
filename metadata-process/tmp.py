""" messing with stuff
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import json, os
import time

time1 = time.time()

# get list of mat files
basedir = jr.getBasedir('NREL')
lmats_fname = [fp for fp in os.listdir(basedir) if 'listmats' in fp][0]
lmats_fpath = os.path.join(basedir,lmats_fname)
with open(lmats_fpath,'r') as f:
    list_mats = json.load(f)
    
i = 50000
h_parms = jr.NRELlistmetadata(i,list_mats)

print(time.time() - time1)