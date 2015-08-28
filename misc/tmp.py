"""
investigating KS test for rhoHat
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
import scipy.io as scio

# load samples of rho_hat
fname = 'rho_hats_KS.mat'
struc = scio.loadmat(fname)
rho_hats = struc['rho_hats']
rhos     = struc['rhos']
mu       = struc['mu']
F_all    = np.arange(1,rho_hats.shape[0]+1)/(rho_hats.shape[0]+1.)

# draw small sample of rho_hats
n_s = 500
n_f = 6000
rho = 0.11
thetas = jr.wrappedCauchySample(n_f,n_s,rho,mu)
samp_rhohats = jr.samplePhaseCoherence(thetas,axis=0)[0]
samp_rhohats = np.sort(samp_rhohats)
F_s = np.arange(1,n_s+1)/(n_s+1.)

plt.clf()
plt.plot(samp_rhohats,F_s)

## get the maximum distance between the CDFs
#x_all = np.concatenate((rhos_1,rhos_2))
#x_int = np.sort(np.unique(x_all))
#F1_int = np.interp(x_all,rhos_1,F_n1)
#F2_int = np.interp(x_all,rhos_2,F_n2)
#D = np.abs(F1_int-F2_int).max()
#
## calculate significance levels
#c = np.array([1.07,1.22,1.36,1.48,1.63,1.73,1.95])
#scale = np.sqrt((n1+n2)/float(n1*n2))
#
#print(D)
#print(c*scale)
#print(D < c*scale)
#
#plt.plot(rhos_2,F_n2)
#
#
## IT WORKSSSSSSSSS    WHOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO