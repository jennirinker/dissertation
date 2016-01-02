"""
Fitting polynomial surface to FAST statistics
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
import os
import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt

# define turbine name and run name
#TurbNames = ['WP0.75A08V00','WP1.5A08V03',
#              'WP3.0A02V02','WP5.0A04V00']
#RunNames = ['Peregrine','TestRun','SmallRun']
TurbName = 'WP5.0A04V00'
RunName  = 'Peregrine'

FigNum = 1

# base directory where the stats are stored
BaseStatDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                'processed_data\\proc_stats'

# -----------------------------------------------------------------------------

# get wind parameters for that run
WindParms = jr.RunName2WindParms(RunName)
URefs, Is, Ls, rhos, n_dups = WindParms['URefs'],WindParms['Is'], \
                              WindParms['Ls'],WindParms['rhos'], \
                              WindParms['n_dups']

# load the stats data
DictPath   = os.path.join(BaseStatDir,RunName,TurbName + '_stats.mat')
stats_dict = scio.loadmat(DictPath,squeeze_me=True)
proc_stats = stats_dict['proc_stats']
calc_stats = [s.rstrip() for s in stats_dict['calc_stats']]
fnames     = [s.rstrip() for s in stats_dict['fnames']]
fields     = [s.rstrip() for s in stats_dict['fields']]
n_fields = len(fields)

# statistic and value to fit RSM to
stat = 'max'
#parm, p_i = 'RootMFlp1', [4,4,4,4]
#parm, p_i = 'RootMOoP1', [4,4,4,4]
#parm, p_i = 'TwrBsMxt', [4,4,4,4]
#parm, p_i = 'TwrBsMyt', [4,4,4,4]
parm, p_i = 'RotTorq', [4,4,2,2]

# extract data for fitting polynomial surface
y = proc_stats[:,calc_stats.index(stat)*n_fields + \
                 fields.index(parm)]                # output data
x = np.empty((y.size,4))
for i_f in range(y.size):
    file_id =  fnames[i_f].rstrip('.out').split('_')[1]
    x[i_f,0] = URefs[int(file_id[0],16)]
    x[i_f,1] = Is[int(file_id[1],16)]
    x[i_f,2] = Ls[int(file_id[2],16)]
    x[i_f,3] = rhos[int(file_id[3],16)]

# fit polynomial surface to all data
coeffs, ps = jr.polyregression(x,y,p_i)  

# ============== plot data and polynomial surface ======================

# --------------------------- U vs. Ti -------------------------------------

i1,i2 = 1,1

# initialize figure
fig = plt.figure(FigNum,figsize=(10,10))
plt.clf()
ax = fig.add_subplot(111, projection='3d')

# mask data to single value of L/rho
mask = np.logical_and(x[:,2] == Ls[i1],x[:,3] == rhos[i2])

# create scatterplot
xplot = x[mask,0]
yplot = x[mask,1]
zplot = y[mask]
ax.scatter(xplot, yplot, zplot, s=9,
           c='k', 
           edgecolors='none')  

# polyfit
X1,X2,X3,X4  = np.meshgrid(URefs,Is,Ls[i1],rhos[i2])
xplot = np.hstack((X1.reshape((X1.size,1)),
                   X2.reshape((X2.size,1)),
                   X3.reshape((X3.size,1)),
                   X4.reshape((X4.size,1))))
A     = jr.myvander(xplot,ps)
z_RSM = np.dot(A,coeffs)
Z     = z_RSM.reshape(X1.shape)
X1 = np.squeeze(X1)
X2 = np.squeeze(X2)
Z = np.squeeze(Z)
ax.plot_wireframe(X1,X2,Z)

# --------------------------- rho vs. L -------------------------------------

i1, i2 = 7,0

# initialize figure
fig = plt.figure(FigNum+1,figsize=(10,10))
plt.clf()
ax = fig.add_subplot(111, projection='3d')

# mask data to single value of L/rho
mask = np.logical_and(x[:,0] == URefs[i1],x[:,1] == Is[i2])

# create scatterplot
xplot = np.log10(x[mask,2])
yplot = x[mask,3]
zplot = y[mask]
ax.scatter(xplot, yplot, zplot, s=9,
           c='k', 
           edgecolors='none')  

# polyfit
X1,X2,X3,X4  = np.meshgrid(URefs[i1],Is[i2],Ls,rhos)
xplot = np.hstack((X1.reshape((X1.size,1)),
                   X2.reshape((X2.size,1)),
                   X3.reshape((X3.size,1)),
                   X4.reshape((X4.size,1))))
A     = jr.myvander(xplot,ps)
z_RSM = np.dot(A,coeffs)
Z     = z_RSM.reshape(X1.shape)
X1 = np.log10(np.squeeze(X3))
X2 = np.squeeze(X4)
Z = np.squeeze(Z)
ax.plot_wireframe(X1,X2,Z)
