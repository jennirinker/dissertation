""" just dicking around
"""
import numpy as np

n_osc = 3
A = 0.6
U = 10.2
T = 630.
f = 1.4
T_steady = 30.
dt = 0.05

t = np.arange(0,T+dt,dt)
u = U*np.ones(t.shape)
i_steady = int(T_steady/dt)
if n_osc is None:
    u[i_steady:] = A*np.sin(2*np.pi*f*(t[i_steady:]-T_steady))
else:
    i_end = i_steady + int(n_osc*(1./f)/dt)
    u[i_steady:i_end] = U+A*np.sin(2*np.pi*f*(t[i_steady:i_end]-T_steady))
    
plt.figure(11)
plt.clf()
plt.plot(t,u)
plt.xlim([29,40])