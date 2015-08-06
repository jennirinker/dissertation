
#import sys
#libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
#if (libpath not in sys.path): sys.path.append(libpath)
#
#import JR_Library.main as jr
#
#basedir = jr.getBasedir('NREL')
#
#list_mats = jr.list_matfiles(basedir,save=1)

import numpy as np
from joblib import Parallel, delayed
import sys

def fcn(i,n):
    return np.array([np.nan,i,n*i])
    
if (__name__ == '__main__'):
    N = 5
    n = 1
#    out = [fcn(i,n) for i in range(N)]
    out = Parallel(n_jobs=2,verbose=1)(delayed(fcn)(i,n) for i in range(N))
    print(out)
