#! /usr/bin/env python
"""
Command-line executable for writing TurbSim input files for Monsoon simulations
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import os

# get turbine name and parameter string from command line
TurbName = sys.argv[1]
FileID   = sys.argv[2]
ParmStr  = sys.argv[3].rstrip('\"').lstrip('\"')


TmplDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\' + \
            'nwtc_python_tools\\templates'
BaseWindDir = '\\\\monsoon-data\\Public\\JRinker\\fast_simulations\\WindDir'

# create TurbSim input dictionary
URef,sig_u,L_u,rho_u,R1,R2 = [float(x) for x in ParmStr.split()]
TSDict = jr.MakeTSDict(TurbName,URef,sig_u,L_u,rho_u,R1,R2)

# Write TurbSim input file
InpName = '{:s}_{:s}.inp'.format(TurbName,FileID)
SpcName = '{:s}_{:s}.spc'.format(TurbName,FileID)
TSDict['SpcName'] = SpcName
WrDir = os.path.join(BaseWindDir,TurbName)
jr.WriteTurbSimInputs(InpName,TSDict,TmplDir,WrDir)