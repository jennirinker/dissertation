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
import statsmodels.api as sm


def StatsErr(x,y,p_i):
    """ Sum-of-squared-error for data x and y with polynomial orders p_i
    """
    
    # perform OLS for all data
    ps_all = jr.GetAllPowers(p_i)
    Xv_all = jr.myvander(x,ps_all)
    results = sm.OLS(y, Xv_all).fit()
    cs_all  = results.params
    
    # extract significant coefficients
    ps_red = ps_all[results.pvalues <= 0.5]
    Xv_red = jr.myvander(x,ps_red)
    
    # fit OLS with reduced coefficients
    cs_red = sm.OLS(y, Xv_red).fit().params
    
    # eliminate insignificant terms
#    cs_all = results.params
#    ps_red = ps[results.pvalues < alpha]
#    cs_red = cs_all[results.pvalues < alpha]
#    Xv_red = jr.myvander(x,ps_red)
#    
#    print(len(cs_all),len(cs_red))
    
    # get reduced model
    yhat   = np.dot(Xv_red,cs_red)
    
    # get residuals and sum
    e    = y - yhat
    err  = np.mean(e ** 2)
    
    return err

# define turbine name and run name
#TurbNames = ['WP0.75A08V00','WP1.5A08V03',
#              'WP3.0A02V02','WP5.0A04V00']
#RunNames = ['Peregrine','TestRun','SmallRun','BigRun2']
TurbName = 'WP5.0A04V00'
RunName  = 'BigRun2'

FigNum = 1

# statistic and value to fit RSM to
stat,parm,units,scale = 'max','RootMFlp1','MN-m',1000.
#stat,parm,units,scale = 'DEL-h','RootMFlp1','MN-m',1000.
#stat,parm,units,scale = 'max','HSShftTq','kN-m',1
#stat,parm,units,scale = 'DEL-h','HSShftTq','kN-m',1
#stat,parm,units,scale = 'max','TwrBsMyt','MN-m',1000
#stat,parm,units,scale = 'DEL-h','TwrBsMyt','MN-m',1000

# base directory where the stats are stored
BaseStatDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                'processed_data\\proc_stats'

alpha = 0.05

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


# extract data for fitting polynomial surface
y = proc_stats[:,calc_stats.index(stat)*n_fields + \
                 fields.index(parm)]                # output data
x = np.empty((y.size,4))
for i_f in range(y.size):
    file_id =  fnames[i_f].rstrip('.out').split('_')[1]
    for i_p in range(len(WindParmsList)):
        x[i_f,i_p] = WindParmsList[i_p][int(file_id[i_p],16)]   # hex to int

# ================= scale input data =====================    

#x = x / np.mean(x,axis=0)

# ================= optimize polynomial coefficients =====================

## parameterize function for data
#ErrFunc = lambda p: StatsErr(x,y,p)
#
#p0 = np.zeros(x.shape[1])
#results = jr.DiscreteOpt(ErrFunc,p0,
#                         verbose=1)
#p_i = results['p_out']
#print(p_i)

#p_i = [7,4,4,4]        # values used in debugging script for Dr. Gavin results
p_i = [7,3,2,2]         # new optimizer output

# ============== plot data and polynomial surface ======================

# get optimal polynomial surface
ps_all = jr.GetAllPowers(p_i)
Xv_all = jr.myvander(x,ps_all)

results = sm.OLS(y, Xv_all).fit()
cs_all  = results.params


pvalues = results.pvalues

ps_red = ps_all[results.pvalues <= 0.5]
Xv_red = jr.myvander(x,ps_red)
cs_red = sm.OLS(y, Xv_red).fit().params

#cs_plot = cs_all
#ps_plot = ps_all
cs_plot = cs_red
ps_plot = ps_red


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
A     = jr.myvander(xplot,ps_plot)
z_RSM = np.dot(A,cs_plot)
Z     = z_RSM.reshape(X1.shape)
X1 = np.squeeze(X1)
X2 = np.squeeze(X2)
Z = np.squeeze(Z)
ax.plot_wireframe(X1,X2,Z)
ax.set_title('{:s} {:s}'.format(stat,parm),fontsize='x-large')
ax.set_xlabel('U',fontsize='x-large')
ax.set_ylabel('I',fontsize='x-large')

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
ax.set_title('{:s} {:s}'.format(stat,parm),fontsize='x-large')
ax.set_xlabel('log10 L',fontsize='x-large')
ax.set_ylabel('rho',fontsize='x-large')

# polyfit ************************************* different than above
X1,X2,X3,X4  = np.meshgrid(WindParmsList[im1][im1i],
                           WindParmsList[im2][im2i],
                           WindParmsList[ip1],
                           WindParmsList[ip2])
xplot = np.hstack((X1.reshape((X1.size,1)),
                   X2.reshape((X2.size,1)),
                   X3.reshape((X3.size,1)),
                   X4.reshape((X4.size,1))))
A     = jr.myvander(xplot,ps_plot)
z_RSM = np.dot(A,cs_plot)
Z     = z_RSM.reshape(X3.shape)
X1 = np.squeeze(X3)                 # ************************
X2 = np.squeeze(X4)                 # ************************
Z = np.squeeze(Z)
ax.plot_wireframe(X1,X2,Z)
