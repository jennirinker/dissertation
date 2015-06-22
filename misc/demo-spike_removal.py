"""
Debug routine for spike removal
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import matplotlib.pyplot as plt

# =============================================================================

# specify data/spike removal parameters to plot
time_tup   = (2012, 2, 13, 17, 0)
ht         = 76
b          = 10.

# load raw data
time_flt   = jr.timetup2flt(time_tup)
t, u, v, w = jr.loadtimeseries('NREL',time_flt,ht,clean=0)

# manually add some spikes for fun
u[:5] = u[6]+2                                          # beginning
u[4000:4005] = u[4000] - 2.5                            # middle down
u[6000:6003] = u[6000]+1.5                              # middle up
u[10700:10701] = u[10700]+0.3                           # middle small up

# clean spike-filled record
u_cl, n_spikes = jr.remove_spikes(u,beta=b)

# plot results
plt.figure(1,figsize=(10,5))
plt.clf()
plt.subplot(2,1,1)
plt.plot(t,u,'r')
plt.plot(t,u_cl)
plt.xlim([-5,605])
plt.ylabel('Wind Velocity')
plt.title(r'Number of spikes: 5, Number removed: {} ($\beta$ = {})' \
    .format(n_spikes,b))

plt.subplot(2,1,2)
plt.plot(t[:-1],u[1:]-u[:-1],'r')
plt.plot(t[:-1],u_cl[1:]-u_cl[:-1])
plt.xlim([-5,605])
plt.ylabel('Wind Velocity Differences')
plt.xlabel('Time')

plt.tight_layout()