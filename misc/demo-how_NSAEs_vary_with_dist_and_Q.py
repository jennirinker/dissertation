"""
Demo of how NSAEs vary with distribution type, quantile.
Conclusion:
    - Not a lot of variation over Q. In general, best dist is smallest for all
      quantiles.
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
import JR_Library.main as jr

import numpy as np
import scipy.io as scio

# %%============================= load data ===================================

# path to matlab-processed metadata table
matpath = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
    '2015-02-28_temporal coherence in data\\code\\tempStruc.mat'

# load the structure
struc = scio.loadmat(matpath)

# extract the information
pdfTable = struc['tempStruc'][0,0][3]

iH = 0
iP = 0

# extract data
x = pdfTable[:,5*iH + iP]
x = np.sort(x)

# %%========================= fit distributions ===============================

# probability distribution candidates
half_cands = ['lognorm','genextreme',\
'chi2','gengamma','exponweib','expon']      # U, sigma_u, L
fine_cands = ['anglit','beta', \
    'genextreme','gengamma','lognorm']      # rho

# quantile threshold values to search
Q_Ts = [0.80,0.85,0.90,0.95,1.00]

# set distribution candidates
if (iP < 3):
    dist_cands = half_cands
elif (iP == 3):
    dist_cands = fine_cands

# initialize error matrix
NSAEs = np.empty((len(dist_cands),len(Q_Ts)))

# initialize parameter list
d_parms = []

# loop through distribution candidates
for iD in range(len(dist_cands)):
    
    # isolate distribution
    dist_name = dist_cands[iD]
    print('Processing {}'.format(dist_name))

    # initialize parameter list
    Q_parms = []

    # loop through threshold values
    for iQ in range(len(Q_Ts)):

        # isolate quantile
        Q_T  = Q_Ts[iQ]
        print('  ...quantile {}'.format(Q_T))

        # calculate threshold value
        x, N = np.sort(x), x.size
        if (Q_T == 1.0): x_T = float('inf')
        else:            x_T  = x[N*Q_T]

        # optimize distribution parameters
        p_main, p_GP = jr.fitcompositeparameters( \
            x,dist_name,x_T)

        # calculate corresponding NSAE
        NSAE = jr.compositeNSAE(x,dist_name,p_main, \
                             x_T,p_GP)

        # save parameters and NSAE
        fit_parms    = (dist_name,p_main,x_T,p_GP,NSAE)
        NSAEs[iD,iQ] = NSAE

        Q_parms.append(fit_parms)

    d_parms.append(Q_parms)

# %%========================= print results ===============================

print('Parm    {:.2f}    {:.2f}    {:.2f}    {:.2f}    {:.2f}'.format(0.8,0.85,0.9,0.95,1))
print('-------------------------------------------------------------------')
for iD in range(len(dist_cands)):
    print('{:<6}  {:.4f}  {:.4f}  {:.4f}  {:.4f}  {:.4f}'.format( \
        dist_cands[iD][:6],
        NSAEs[iD,0],NSAEs[iD,1],NSAEs[iD,2],NSAEs[iD,3],NSAEs[iD,4]))



