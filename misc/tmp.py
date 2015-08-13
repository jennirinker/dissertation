"""
testing routine to loop through M4 files, try to load them,
save a list of unloadable files, then re-download those files and try again
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr

# gat base directory for data
basedir = jr.getBasedir('NREL')
