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
#RunNames = ['Peregrine','TestRun','SmallRun','BigRun2']
TurbName = 'WP5.0A04V00'
RunName  = 'BigRun2'

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
    

# ================= optimize polynomial coefficients =====================

## parameterize function for data
#ErrFunc = lambda p: StatsErr(x,y,p)
#
#p0 = np.zeros(x.shape[1])
#results = jr.DiscreteOpt(ErrFunc,p0,
#                         verbose=1)
#p_i = results['p_out']
#print(p_i)

# ============== plot residuals ======================

# get optimal polynomial surface
betas, ps = jr.polyregression(x,y,p_i)
A         = jr.myvander(x,ps)
yhat      = np.dot(A,betas)
maxerr    = np.abs(yhat - y).max()


n_x, n_y = len(WindParmsList[2]),len(WindParmsList[3])
xplot = (np.arange(n_x)/float(n_x))*0.9 + 0.1
yplot = (np.arange(n_y)/float(n_y))[::-1]*0.9 + 0.08
dx,dy = 0.75*1./n_x, 0.75*1./n_y

plt.figure(1,figsize=(10,10))
plt.clf()

# loop through Ls and rhos
for i_x in range(n_x):
    for i_y in range(n_y):
        
        # extract data
        mask = np.logical_and(x[:,2] == WindParmsList[2][i_x],
                              x[:,3] == WindParmsList[3][i_y])
        x_mask = x[mask,:]
        
        # get average value of data
        x1, x2 = np.unique(x_mask[:,0]), np.unique(x_mask[:,1])
        n1, n2 = len(x1), len(x2)
        n_uniq = n1*n2
        x_avg, y_avg = np.empty((n_uniq,x.shape[0])), np.empty(n_uniq)
        for i1 in range(n1):
            for i2 in range(n2):
                
                
        
        
        plt.axes([xplot[i_x],yplot[i_y],dx,dy])



