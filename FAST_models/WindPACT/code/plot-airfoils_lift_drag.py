"""
Plot lift and drag coefficients of specified airfoil
"""
import numpy as np

fname = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
            'dissertation\\FAST_models\\FAST7\\WP1.5A08V03\\AeroData\\' + \
            's818_2703.dat'

AFData = np.empty((1e4,3))
with open(fname,'r') as f:
    i_line = 0
    i_AF   = 0
    for line in f:
        data = line.rstrip('\n').strip()
        if data and (i_line > 14):
            AFData[i_AF,:] = [float(l) for l in line.split()[:3]]
            i_AF += 1
        i_line += 1

AFData = AFData[:i_AF]

plt.figure(1)
plt.clf()
plt.plot(AFData[:,0],AFData[:,1],label='C_L_818')
plt.plot(AFData[:,0],AFData[:,2],label='C_D_818')

fname = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
            'dissertation\\FAST_models\\FAST7\\WP0750_v7\\AeroData\\' + \
            's825_2102.dat'

AFData = np.empty((1e4,3))
with open(fname,'r') as f:
    i_line = 0
    i_AF   = 0
    for line in f:
        if i_line > 11:
            AFData[i_AF,:] = [float(l) for l in line.split()[:3]]
            i_AF += 1
        i_line += 1

AFData = AFData[:i_AF]


plt.plot(AFData[:,0],AFData[:,1],label='C_L_825')
plt.plot(AFData[:,0],AFData[:,2],label='C_D_825')
plt.legend()


fname = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
            'dissertation\\FAST_models\\FAST7\\WP0750_v7\\AeroData\\' + \
            's826_1602.dat'

AFData = np.empty((1e4,3))
with open(fname,'r') as f:
    i_line = 0
    i_AF   = 0
    for line in f:
        if i_line > 11:
            AFData[i_AF,:] = [float(l) for l in line.split()[:3]]
            i_AF += 1
        i_line += 1

AFData = AFData[:i_AF]


plt.plot(AFData[:,0],AFData[:,1],label='C_L_826')
plt.plot(AFData[:,0],AFData[:,2],label='C_D_826')
plt.legend()