"""
Summarize processed wind parameters from a given dataset
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr

# choose which dataset
dataset = 'NREL'            # NREL, NREL-mat, fluela
