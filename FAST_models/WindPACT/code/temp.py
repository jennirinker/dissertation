"""
Temp testing reading text files copied from excel
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr



fname = 'textfiles//0.75_Blade.txt'
header = 1
units  = 0

data, header = jr.mygenfromtxt(fname)