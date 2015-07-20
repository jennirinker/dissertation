"""
Debug routine for spike removal
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import matplotlib.pyplot as plt
import random

# %%===========================================================================

# specify data/spike removal parameters to plot
time_tup   = (2012, 2, 13, 17, 0)
ht         = 76
b          = 10.

# load raw data
time_flt   = jr.timetup2flt(time_tup)
outdict = jr.loadtimeseries('NREL','Sonic_u_15m',time_flt)
t, x, x_cl = outdict['time'], outdict['raw'], outdict['clean']

# %%===========================================================================
# testing spike removal algorithm

# manually add some spikes for fun
x[:5] = x[6]+2                                          # beginning
x[4000:4005] = x[4000] - 2.5                            # middle down
x[6000:6003] = x[6000]+1.5                              # middle up
x[10700:10701] = x[10700]+0.4                           # middle small up

# clean spike-filled record
u_cl = jr.remove_spikes(x,beta=b)

# plot results
plt.figure(1,figsize=(10,5))
plt.clf()

# wind velocity
plt.subplot(2,1,1)
plt.plot(t,x,'r')
plt.plot(t,u_cl)
plt.xlim([-5,605])
plt.ylabel('Wind Velocity')
plt.title(r'Number of spikes: 5 ($\beta$ = {})' \
    .format(b))

# wind velocity differences
plt.subplot(2,1,2)
plt.plot(t[:-1],x[1:]-x[:-1],'r')
plt.plot(t[:-1],u_cl[1:]-u_cl[:-1])
plt.xlim([-5,605])
plt.ylabel('Wind Velocity Differences')
plt.xlabel('Time')

plt.tight_layout()

# %%===========================================================================
# testing quantization detection

# manually add quatnization occurrences
n_quants = 11
n_zeros  = 4
i_quant = random.sample(range(1,6000,n_zeros),n_quants)
#i_quant = [3297, 5725, 201, 4389, 3973, 2433, 3605, 3085, 3801, 5433, 797]
x_quant = np.copy(x_cl)
for i in range(n_quants):
    x_quant[i_quant[i]:i_quant[i]+n_zeros-1+1] = x_quant[i_quant[i]-1]

flag = jr.is_quantized(x_quant)
print('quant locs: ',sorted(i_quant))
print(flag)



