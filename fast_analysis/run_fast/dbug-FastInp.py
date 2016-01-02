#! /usr/bin/env python
"""
Command-line executable for writing FAST input files for Monsoon simulations
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\nwtc_python_tools'
if (libpath not in sys.path): sys.path.append(libpath)

import os, jr_fast

# get turbine name and parameter string from command line
#TurbName    = sys.argv[1]
#WindPath = sys.argv[2]

TurbName = 'WP0.75A08V00'
WindPath = '\\\\monsoon-data\\Public\\JRinker\\fast_simulations\\' + \
            'WindDir\\WP0.75A08V00\\WP0.75A08V00_00000.bts'
            
# output and input directories
FastDir = os.path.join('\\\\monsoon-data\\Public\\JRinker\\' + \
                       'fast_simulations\\FastDir',TurbName)
ModlDir = os.path.join('\\\\monsoon-data\\Public\\JRinker\\' + \
                       'fast_simulations\\ModlDir',TurbName)


# get name for .fst file
FastName = os.path.splitext(os.path.split(WindPath)[-1])[0]

# write FAST file
jr_fast.WriteFastADOne(TurbName,WindPath,FastName,ModlDir,FastDir)
