"""
From a list of distribution candidates, fit single distribution to data and
save top four distributions.
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
import JR_Library.main as jr

import scipy.io as scio
import numpy as np
import matplotlib.pyplot

# %%============================= load data ===================================

# path to matlab-processed metadata table
matpath = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
    '2015-02-28_temporal coherence in data\\code\\tempStruc.mat'

# load the structure
struc = scio.loadmat(matpath)

# extract the information
pdfTable = struc['tempStruc'][0,0][3]

# %%========================= fit distributions ===============================

# probability distribution candidates
half_cands = ['lognorm','genextreme',\
    'chi2','gengamma','exponweib','expon']      # half-infinite support
fine_cands = ['anglit','arcsine','beta','cosine', \
    'genextreme','vonmises','wrapcauchy']       # finite support
x_T = float('inf')                              # fitting single distributions
                    
iP = 0
iH = 0

# extract data
x = pdfTable[:,5*iH + iP]
x = np.sort(x)

# choose half-infinite or finite support distributions
if (iP < 3): dist_cands = half_cands
elif (iP ==3): dist_cands = fine_cands

d_parms = []

# loop through distribution candidates
for iD in range(len(dist_cands)):
        
    # separate distribution name
    dist_name = dist_cands[iD]
    print('    Processing {}'.format(dist_name))
    #    
    # fit distribution to data
    p_main, p_GP = jr.fitcompositeparameters(x,dist_name,x_T)

    # calculate NSAE
    NSAE = jr.compositeNSAE(x,dist_name,p_main,x_T,p_GP)
    
    # create parameter list
    fit_parms = [dist_name,p_main,x_T,p_GP,NSAE]
    
    d_parms.append(fit_parms)

# %% ======================= plot distributions ===============================

plt.figure(1)
plt.clf()

N = x.size
F_emp = np.arange(N)/(N+1.)
plt.plot(x,F_emp)

legstr = ['empirical (0.0)']
for iD in range(len(dist_cands)):
    F = jr.compositeCDF(x,d_parms[iD][0],d_parms[iD][1], \
        d_parms[iD][2],d_parms[iD][3])
    plt.plot(x,F)

    leg_entry = d_parms[iD][0] + ', ({:.4f})'.format(d_parms[iD][4])
    legstr = legstr + [leg_entry]

plt.legend(legstr,loc='lower right')



