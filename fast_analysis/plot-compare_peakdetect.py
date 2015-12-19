"""
Comparison of peak detection algorithms for WISDEM algorithm versus algorithm
found online at gist github
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
import matplotlib.pyplot as plt
import numpy as np
from WISDEM_rainflow.rainflow import determine_peaks
from JR_Library.peakdetect import peakdetect

# make plots pretty
plt.style.use('C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
                '2016-02-15_dissertation\\figure_code\\' + \
                'duke_presentation.mplstyle')

# choose whether to reload FAST dictionary
reload_dict = 1

# time to look ahead for gist algorithm
t_lookahead = 0.5

# what value to analyze
key = 'RootMOoP1'

# path to FAST file
fpath_fast = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation\\' + \
                'FAST_models\\FAST7\\WP0.75A08V00_equil\\WP0.75A08V00_24134.out'

# reload dictionary if requested
if reload_dict:
    fast_dict = jr.ReadFASTFile(fpath_fast)
    
# get time series from dictionary
data = fast_dict['Data']
fields = fast_dict['Fields']
units = fast_dict['Units']
t = data[:,fields.index('Time')]
y = data[:,fields.index(key)]

# number of steps to look ahead for Gist algorithm
lookahead = int(t_lookahead / t[1])

# get peak values and indices using WISDEM algorithm
peaks, ipeaks = determine_peaks(y)
xpeaks_W = t[ipeaks]
ypeaks_W = y[ipeaks]

# get peak values and indices using gist algorithm
maxinfo, mininfo = peakdetect(y, t, 
                              lookahead=1)
xpeaks_g = np.array([tup[0] for tup in maxinfo] + [tup[0] for tup in mininfo])
ypeaks_g = np.array([tup[1] for tup in maxinfo] + [tup[1] for tup in mininfo])

# resort gist algorithm into increasing x value
isort    = xpeaks_g.argsort()
xpeaks_g = xpeaks_g[isort]
ypeaks_g = ypeaks_g[isort]

# -------------------------- plot results -------------------------------------
plt.figure(1,figsize=(6.5,5.))
plt.clf()

ms = 8
ax = plt.subplot(211)
plt.plot(t,y,'k',
         lw=2)
plt.plot(xpeaks_W, ypeaks_W, '^',
         ms=ms,mfc='none',mec='r',mew=2,label='WISDEM')
plt.plot(xpeaks_g, ypeaks_g, 'o',
         ms=ms,mfc='none',mec='b',mew=2,label='Gist - Github')

plt.xlim([72,77])
plt.ylim([200,700])
plt.xlabel('Time (s)')

ax.legend(loc='lower right', bbox_to_anchor=(1.0, 1.05))        
         
plt.subplot(212)
plt.plot(t,y,'k',
         lw=2)
plt.plot(xpeaks_W, ypeaks_W,'r',
         lw=2)
plt.plot(xpeaks_g, ypeaks_g,'b',
         lw=2)

plt.xlim([72,77])
plt.ylim([200,700])
plt.xlabel('Time (s)')

plt.figtext(0.05,0.5,key + ' ' + units[fields.index(key)],
            ha='center',va='center',rotation='vertical',
            fontsize='x-large')

plt.tight_layout(rect=[0.07,0,1.,0.91])


