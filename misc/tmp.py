# -*- coding: utf-8 -*-
"""
Created on Wed Jul 01 07:25:38 2015

@author: jrinker
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import scipy.io as scio

dataset = 'NREL'

idx_high = metadata[:,8].argmax()
time_flt = metadata[idx_high,0]
time_tup = jr.timeflt2tup(time_flt)
ht = metadata[idx_high,2]
#print(jr.loadtimeseries(dataset,'Cup_WS_10m',time_hlt)['flags'])
fpath = jr.NRELtime2fpath(time_flt)
struc20 = scio.loadmat(fpath)
jr.calculatefield(dataset,struc20,ht)
t = np.arange(12000)*0.05


plt.figure(1)
plt.clf()
plt.subplot(411)
datfield = jr.field2datfield(dataset,'Sonic_u',ht)
out = jr.loadtimeseries(dataset,datfield,time_flt)
plt.plot(t,out['raw'])
#plt.plot(out['clean'])
plt.title('Sonic_u,  (y,m,d,h,min,ht) = {}'.format(time_tup+(ht,)))

plt.subplot(412)
plt.plot(t,out['raw'])
#plt.plot(out['clean'])
plt.xlim([400,600])
plt.title('Sonic_u (zoomed)'.format(time_tup+(ht,)))

plt.subplot(413)
datfield = jr.field2datfield(dataset,'Precipitation',ht)
out = jr.loadtimeseries('NREL',datfield,time_flt)
plt.plot(t,out['raw'])
#plt.plot(out['clean'])
plt.title('Precipitation')

plt.subplot(414)
datfield = jr.field2datfield(dataset,'Wind_Speed_Cup',88)
out = jr.loadtimeseries(dataset,datfield,time_flt)
plt.plot(t,out['raw'])
#plt.plot(out['clean'])
plt.title('Cup Wind Speed')

plt.tight_layout()