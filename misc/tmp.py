"""
investigating KS test for rhoHat
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np

# calculate first CDF
#n1  = 10000
#n_f = 6000
#rho = 0.11
#mu  = np.pi
#thetas1 = jr.wrappedCauchySample(n_f,n1,rho,mu)
#rho_hats1 = jr.samplePhaseCoherence(thetas1,axis=0)[0]

#F_n1 = np.arange(1,n1+1)/(n1+1.)
#rhos_1 = np.sort(rho_hats1)

# calculate second CDF
n2  = 500
n_f = 6000
rho = 0.11
mu  = np.pi
thetas2 = jr.wrappedCauchySample(n_f,n2,rho,mu)
rho_hats2 = jr.samplePhaseCoherence(thetas2,axis=0)[0]

F_n2 = np.arange(1,n2+1)/(n2+1.)
rhos_2 = np.sort(rho_hats2)

# get the maximum distance between the CDFs
x_all = np.concatenate((rhos_1,rhos_2))
x_int = np.sort(np.unique(x_all))
F1_int = np.interp(x_all,rhos_1,F_n1)
F2_int = np.interp(x_all,rhos_2,F_n2)
D = np.abs(F1_int-F2_int).max()

# calculate significance levels
c = np.array([1.07,1.22,1.36,1.48,1.63,1.73,1.95])
scale = np.sqrt((n1+n2)/float(n1*n2))

print(D)
print(c*scale)
print(D < c*scale)

plt.plot(rhos_2,F_n2)


# IT WORKSSSSSSSSS    WHOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO