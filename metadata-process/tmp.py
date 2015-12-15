"""
testing reading texas tech data with .mat and text files
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import scipy.signal
import numpy as np
import matplotlib.pyplot as plt
from JR_Library.peakdetect import peakdetect

i = 10000
x = np.linspace(0,3.7*np.pi,i)
y = (0.3*np.sin(x) + np.sin(1.3 * x) + 0.9 * np.sin(4.2 * x) + 0.00 * \
        np.random.randn(i))
y *= -1

_max, _min = peakdetect(y, x,
                        lookahead=300,delta=0.30)
xmax = [p[0] for p in _max]
ymax = [p[1] for p in _max]
xmin = [p[0] for p in _min]
ymin = [p[1] for p in _min]

plt.figure(1)
plt.clf()

plt.plot(x,y)
plt.plot(xmax,ymax,'ko')
plt.plot(xmin,ymin,'ko')
#plt.ylim([-1.5,1.5])


