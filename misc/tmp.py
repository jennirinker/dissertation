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

time_flt = metadata[3,0]
print(jr.loadtimeseries(dataset,'Cup_WS_10m',time_hlt)['flags'])
fpath = jr.NRELtime2fpath(time_flt)
struc20 = scio.loadmat(fpath)
ht = 15
jr.calculatefield(dataset,struc20,ht)


plt.figure(1)
plt.clf()
plt.subplot(211)
datfield = jr.field2datfield(dataset,'Sonic_T',ht)
out = jr.loadtimeseries(dataset,datfield,time_flt)
plt.plot(out['raw'])
plt.plot(out['clean'])

plt.subplot(212)
datfield = jr.field2datfield(dataset,'Precipitation',ht)
out = jr.loadtimeseries('NREL',datfield,time_flt)
plt.plot(out['raw'])
plt.plot(out['clean'])