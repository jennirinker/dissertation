"""
Testing discrete optimization for determining 
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
import numpy as np

def CalcError(p):
    p_ex = np.array([2,9.3,3])
    err = np.sum((p - p_ex)**2)
    return err
    

# iteration parameters
p0 = np.array([1,1,1])

results = jr.DiscreteOpt(CalcError,p0,
                max_iters=100,err_thresh=1e-16,
                derr_thresh=0,n_last=4,
                verbose=0)
            
print('Number of iterations: {:d}'.format(results['num_iters']))
print('Output parameters: [{:s}]'.format(','.join([str(x) for x in results['p_out']])))
print('Error: {:.2f}'.format(results['error']))