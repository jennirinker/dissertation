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
from mpl_toolkits.mplot3d import Axes3D


def StatsErr(x,y,p_i):
    """ Sum-of-squared-error for data x and y with polynomial orders p_i
    """
    
    # perform linear regression
    betas, ps = jr.polyregression(x,y,p_i) 
    
    # calculate model values
    X     = jr.myvander(x,ps)
    yhat  = np.dot(X,betas)
    
    # get residuals and sum
    e    = y - yhat
#    err  = np.mean(e ** 2)
    err  = np.mean(np.abs(e))
    
    return err

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
WindParmsList = [URefs,Is,np.log10(Ls),rhos]

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
parm = 'RootMFlp1'
#parm, p_i = 'RootMOoP1', [4,4,4,4]
#parm, p_i = 'TwrBsMxt', [4,4,4,4]
#parm, p_i = 'TwrBsMyt', [4,4,4,4]
#parm, p_i = 'RotTorq', [4,4,2,2]

# extract data for fitting polynomial surface
y = proc_stats[:,calc_stats.index(stat)*n_fields + \
                 fields.index(parm)]                # output data
x = np.empty((y.size,4))
for i_f in range(y.size):
    file_id =  fnames[i_f].rstrip('.out').split('_')[1]
    for i_p in range(len(WindParmsList)):
        x[i_f,i_p] = WindParmsList[i_p][int(file_id[i_p],16)]   # hex to int
#    x[i_f,1] = Is[int(file_id[1],16)]
##    x[i_f,2] = Ls[int(file_id[2],16)]
#    x[i_f,2] = np.log10(Ls[int(file_id[2],16)])
#    x[i_f,3] = rhos[int(file_id[3],16)]
    

# ================= optimize polynomial coefficients =====================

## parameterize function for data
#ErrFunc = lambda p: StatsErr(x,y,p)
#
#p0 = np.zeros(x.shape[1])
#results = jr.DiscreteOpt(ErrFunc,p0,
#                         verbose=1)
#p_i = results['p_out']
#print(p_i)

# ============== plot data and polynomial surface ======================

# get optimal polynomial surface
betas, ps = jr.polyregression(x,y,p_i)

# --------------------------- U vs. Ti -------------------------------------

# indices for plotting and masking
ip1,ip2   = 0, 1                # plot U and Ti
im1,im2   = 2, 3                # mask by L and rho
im1i,im2i = 1, 1                # mask by 2nd L and rho values
FigNum    = 1                   # figure number

# initialize figure
fig = plt.figure(FigNum,figsize=(10,10))
plt.clf()
ax = fig.add_subplot(111, projection='3d')

# mask data to single value of L/rho
mask = np.logical_and(x[:,im1] == WindParmsList[im1][im1i],
                      x[:,im2] == WindParmsList[im2][im2i])
#mask = np.logical_and(x[:,2] == Ls[i1],x[:,3] == rhos[i2])

# create scatterplot
xplot = x[mask,ip1]
yplot = x[mask,ip2]
zplot = y[mask]
ax.scatter(xplot, yplot, zplot, s=9,
           c='k', 
           edgecolors='none')  

# polyfit
#X1,X2,X3,X4  = np.meshgrid(URefs,Is,Ls[i1],rhos[i2])
X1,X2,X3,X4  = np.meshgrid(WindParmsList[ip1],WindParmsList[ip2],
                           WindParmsList[im1][im1i],
                           WindParmsList[im2][im2i])
xplot = np.hstack((X1.reshape((X1.size,1)),
                   X2.reshape((X2.size,1)),
                   X3.reshape((X3.size,1)),
                   X4.reshape((X4.size,1))))
A     = jr.myvander(xplot,ps)
z_RSM = np.dot(A,betas)
Z     = z_RSM.reshape(X1.shape)
X1 = np.squeeze(X1)
X2 = np.squeeze(X2)
Z = np.squeeze(Z)
ax.plot_wireframe(X1,X2,Z)

# --------------------------- rho vs. L -------------------------------------

# indices for plotting and masking
ip1,ip2   = 2, 3
im1,im2   = 0, 1
im1i,im2i = 3, 0
FigNum    = 2                   # figure number

# initialize figure
fig = plt.figure(FigNum,figsize=(10,10))
plt.clf()
ax = fig.add_subplot(111, projection='3d')

# mask data to single value of L/rho
mask = np.logical_and(x[:,im1] == WindParmsList[im1][im1i],
                      x[:,im2] == WindParmsList[im2][im2i])

# create scatterplot
xplot = x[mask,ip1]
yplot = x[mask,ip2]
zplot = y[mask]
ax.scatter(xplot, yplot, zplot, s=9,
           c='k', 
           edgecolors='none')  

# polyfit ************************************* different than above
X1,X2,X3,X4  = np.meshgrid(WindParmsList[im1][im1i],
                           WindParmsList[im2][im2i],
                           WindParmsList[ip1],
                           WindParmsList[ip2])
xplot = np.hstack((X1.reshape((X1.size,1)),
                   X2.reshape((X2.size,1)),
                   X3.reshape((X3.size,1)),
                   X4.reshape((X4.size,1))))
A     = jr.myvander(xplot,ps)
z_RSM = np.dot(A,betas)
Z     = z_RSM.reshape(X3.shape)
X1 = np.squeeze(X3)                 # ************************
X2 = np.squeeze(X4)                 # ************************
Z = np.squeeze(Z)
ax.plot_wireframe(X1,X2,Z)
