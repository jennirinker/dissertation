"""
Debugging fatigue DEL calculations
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import scipy.signal
import numpy as np
import matplotlib.pyplot as plt
from JR_Library.rainflow.rainflow import determine_peaks
import JR_Library.main as jr
from JR_Library.peakdetect import peakdetect

reload = 0

fpath_fast = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation\\' + \
                'FAST_models\\FAST7\\WP0.75A08V00_equil\\WP0.75A08V00_24134.out'

if reload:

    fast_dict = jr.ReadFASTFile(fpath_fast)
    data = fast_dict['Data']
    fields = fast_dict['Fields']
    units = fast_dict['Units']
    
    x = data[:,fields.index('Time')]
    y = data[:,fields.index('RootMOoP1')]

lookahead = 0.5/x[1]
peaks1, ipeaks1 = determine_peaks(y)
maxs2, mins2 = peakdetect(y,x,
                          lookahead=lookahead)

xpeaks1 = x[ipeaks1]
ypeaks1 = y[ipeaks1]

xpeaks2 = [tup[0] for tup in maxs2] + [tup[0] for tup in mins2]
ypeaks2 = [tup[1] for tup in maxs2] + [tup[1] for tup in mins2]


plt.figure(1)
plt.clf()

plt.plot(x,y)
plt.plot(xpeaks1,ypeaks1,'o',
         ms=10,mfc='none',mec='r',mew=2,label='WISDEM')
plt.plot(xpeaks2,ypeaks2,'o',
         ms=10,mfc='none',mec='g',mew=2,label='Gist')
plt.xlim([70,80])
#plt.plot(xmin,ymin,'ko')
#plt.ylim([-1.5,1.5])


