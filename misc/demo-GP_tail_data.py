"""
Demonstration of fitting GP distribution to tail data using sigma_u @ 15 m.
Accomplishes several things:
    1)  Verifies that <fitcompositedistribution> in JR_Library works for 
        fitting single distributions.
    2)  Shows the NSAE values for single and composite distributions (for this
        example it seems like they are very similar?)
    3)  Plots CDF, CDF errors for single and composite distributions.
Conclusions:
    - Fitting method works well, quickly
    - Simultaneously fitting main and GP parameters works
    - At least for this selection of data, use of GP might not be ideal
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import numpy as np
import JR_Library.main as jr
import scipy.io as scio
import matplotlib.pyplot as plt

def fitsingleparameters(x,dist_name):
    import sys
    libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
    if (libpath not in sys.path): sys.path.append(libpath)
    import JR_Library.main as jr
    import scipy.stats
    from scipy.optimize import minimize
    import numpy as np

    # initialize CDFs
    dist_main = getattr(scipy.stats, dist_name)

    # get initial guess using MLE
    p_main0 = dist_main.fit(x)

    # define error function
    fun = lambda p_main: jr.compositeNSAE(x,dist_name,p_main)

    # perform optimization
    res = minimize(fun,p_main0)

    # convert numpy array to tuple
    param_opt = tuple(np.ndarray.tolist(res.x))

    return param_opt

# ============================== extract data =================================

# path to matlab-processed metadata table
matpath = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
    '2015-02-28_temporal coherence in data\\code\\tempStruc.mat'

# load the structure
struc = scio.loadmat(matpath)

# extract the information
pdfTable = struc['tempStruc'][0,0][3]

# extract data *** sigma
x = pdfTable[:,1]
x = np.sort(x)

# miscellaneous parameters
dist_name = 'lognorm'
N   = x.size
F_emp = np.arange(1,N+1)/(N+1.)

# ======================== fit single distribution ============================

# test function above
p_opt_single = fitsingleparameters(x,dist_name)

# test composite function for fitting single
p_opt_comp = jr.fitcompositeparameters(x,dist_name)[0]

# print parameters
print('Verification that JR_Library function works for single distributions')
print('p_opt_single = {}').format(p_opt_single)
print('p_opt_comp = {}').format(p_opt_comp)

# single CDF
F_single = jr.compositeCDF(x,dist_name,p_opt_comp)

# ====================== fit composite distribution ===========================

# get threshold value for 85%
Q = 0.85
x_T = x[Q*N]

# optimize parameters
p_main_opt, p_GP_opt = jr.fitcompositeparameters(x,dist_name,x_T)

# print composite distribution parameters for fun
print('\nComposite distribution parameters:')
print('Main: {}'.format(p_main_opt))
print('GP: {}'.format(p_GP_opt))

# CDFs
F_comp = jr.compositeCDF(x,dist_name,p_main_opt,x_T,p_GP_opt)

# ========================== compare NSAE values ==============================
NSAE_single = jr.compositeNSAE(x,dist_name,p_opt_comp)
NSAE_comp   = jr.compositeNSAE(x,dist_name,p_main_opt,x_T,p_GP_opt)

print('\nNSAE errors:')
print('Single: {}'.format(NSAE_single))
print('Comp:   {}'.format(NSAE_comp))

# ====================== plot comparison of results ===========================

# make figure
plt.figure(1)
plt.clf()

# subplots of single CDF/CDF error
plt.subplot(2,2,1)
plt.plot(x,F_emp)
plt.plot(x,F_single)
plt.title('Single Distribution')

plt.subplot(2,2,2)
plt.plot(x,F_single-F_emp)
plt.plot(x,F_comp-F_emp)
plt.plot(x[Q*N],F_comp[Q*N]-F_emp[N*Q],'o',mec='r',mfc='none',markersize=8,mew=2)
plt.title('CDF Errors')

# subplots of composite CDF/CDF error
plt.subplot(2,2,3)
plt.plot(x,F_emp)
plt.plot(x,F_comp)
plt.plot(x[Q*N],F_comp[Q*N],'o',mec='r',mfc='none',markersize=8,mew=2)
plt.title('Composite Distribution')

plt.subplot(2,2,4)
plt.plot(x,F_single-F_comp)
plt.plot(x[Q*N],F_single[Q*N]-F_comp[Q*N],'o',mec='r',mfc='none',markersize=8,mew=2)
plt.title('Difference in Fit CDFs')





