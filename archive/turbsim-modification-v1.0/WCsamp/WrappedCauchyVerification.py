""" Calls WrappedCauchySample executable to write sample to file, then compares
written results with theoretical expectations.
"""

import os
import numpy as np
import sys
sys.path.append('../..')
import JR_Library.misc as misc
import matplotlib.pyplot as plt

n = 5000            # number of points in sample
rho = 0.2           # concentration parameter
mu = np.pi/4        # location parameter

# filename where output is stored
fname = "WrappedCauchyOut.txt"

# call Fortran executable
os.system("./WrappedCauchySample {} {} {}". \
          format(n,rho,mu))

# load Wrapped Cauchy output into numpy array
WC_out = misc.wrap(np.loadtxt(fname))

# theory
thetaPlot = np.linspace(-np.pi,np.pi,100)
f = misc.wrappedCauchyPDF(thetaPlot,rho,mu)

# create figure and axes
plt.figure(1,figsize=(10,8))
plt.clf()
plt.axes([0.10,0.10,0.85,0.82])

# plot histogram and smooth line
nbins = 20
plt.hist(WC_out,bins=nbins,normed=True, \
         histtype='step',linewidth=2.,\
         label='Fortran sample')# plot results
plt.plot(thetaPlot,f,linewidth=2.,label='Theory')
plt.xlabel('Angle (rad)')
plt.ylabel('Frequency')
plt.legend()
titleStr = r'$N = {:.0f}$, $\rho = {:.2f}$, $\mu = {:.2f}$'.format(n,rho,mu)
plt.title(titleStr)

# show plot
plt.ion()
plt.show()

