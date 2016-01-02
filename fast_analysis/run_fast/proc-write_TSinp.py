#! /usr/bin/env python
"""
Write TurbSim input files
"""
import sys,os
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import random

# define turbine name
TurbNames = ['WP0.75A08V00','WP1.5A08V03','WP3.0A02V02','WP5.0A04V00']
BaseWindDir = '\\\\monsoon-data\\Public\\JRinker\\fast_simulations\\WindDir'


# wind parameters
#URefs  = [5.0, 7.0, 9.0, 10.0, 10.5, 11.0, 11.5, \
#        12.0, 13.0, 16.0, 19.0, 22.0]
#Is     = [0.1, 0.2, 0.3, 0.4, 0.5]
#logLs  = [1.5, 2.0, 2.5, 3.0]
#rhos   = [0.0, 0.1, 0.2, 0.3, 0.4]
#n_dups = 5
URefs  = [9.0, 10.0, 10.5, 11.0, 11.5, 12.0]
Is     = [0.1, 0.2, 0.3, 0.4, 0.5]
logLs  = [2.0]
rhos   = [0.4]
n_dups = 10

# loop through turbines
for TurbName in TurbNames:
    
    print('Creating TurbSim input files for turbine \"{:s}\"'.format(TurbName))

    # define output directory for TurbSim files
    out_dir = os.path.join(BaseWindDir,TurbName)
    
    # create output directory if it doesn't exist
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    
    # write TurbSim file
    TSDict = jr.MakeTSDict(TurbName,URef,sig_u,L_u,rho_u,R1,R2)