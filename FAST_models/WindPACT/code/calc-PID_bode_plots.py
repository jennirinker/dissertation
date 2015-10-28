"""
Figure out natural frequencies, damping, of each section
"""
import numpy as np


a = 0.81
tau = 7.917**(-1)

wns = 2*np.pi*np.logspace(-2,2,100)
#C = (1/a*1/tau + 1j*wns)/(1/tau + 1j*wns)

a0,a1 = 12.32,1.4
b0,b1 = 7.917,1.0

#C = (a0 + a1*(1j*wns))/(b0 + b1*(1j*wns))
#C = 

plt.figure(10)
plt.clf()

plt.subplot(2,1,1)
plt.loglog(wns/2/np.pi,np.abs(C))
plt.axis('tight')

plt.subplot(2,1,2)
plt.semilogx(wns/2/np.pi,np.angle(C))