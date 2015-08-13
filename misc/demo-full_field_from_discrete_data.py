"""
Constructing a full turbulent field from a discrete set of "measured" time 
series.
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
from scipy import optimize
import scipy.io as scio

# **************** TO DO ****************************
# -- add in check to make sure grid does not overlap points
#     (if it does, just use time series)

# ============= inputs =============

# define locations of data points
y_data = np.array([0.,0.,0.,0.,0.,0.])
z_data = np.array([15.,30,50,76,100,131.])

# define desired grid
y_grid = np.linspace(-15,15,9)
z_grid = np.linspace(20,160,10)

# ============= create or load time series =============
timestamp = (2013,11,3,14,30)
fpath     = 'C:\\Users\\jrinker\\Dropbox\\research\\temp_M4_files\\11_03_2013_14_30_00_026.mat'
datfields = ['Sonic_u_15m','Sonic_u_30m','Sonic_u_50m','Sonic_u_76m',\
                'Sonic_u_100m','Sonic_u_131m']
n_t, dt = 12000, 0.05
t_data = np.arange(n_t)*dt
struc = scio.loadmat(fpath)
n_data = y_data.size
u_data = np.empty((n_t,n_data))
for i in range(len(datfields)):
    x_raw = np.squeeze(struc[datfields[i]][0,0][0])  
    x_cl = jr.cleantimeseries(t,x_raw)
    u_data[:,i] = x_cl
#n_t, dt = 12000, 0.05
#zhub, Vhub, turbc = 90., 8., 'A'
#U_z  = jr.IEC_VelProfile(z_data,zhub,Vhub)
#Iref = jr.IEC_Iref(turbc)
#sig1 = jr.IEC_Sigma1(Iref,Vhub)
#L = 8.1 * 42
#u_data = np.empty((n_t,n_data))
#for i in range(n_data):
#    u_data[:,i] = np.squeeze(jr.generateKaimal1D(n_t,1,dt, \
#        U_z[i],sig1,L/U_z[i],0.0,0.0)[1])
        

    
# ============= spatially correlate phases, create field =============


# %%  ============= test plot =============
#
# create plot
plt.figure(1,figsize=(6,6))
plt.clf()
ax = plt.axes([0.12,0.12,0.81,0.66])
alpha1 = 1
alpha2 = (z_grid[-1]-z_grid[0])/(y_grid[-1]-y_grid[0])

# plot data phases
for i in range(n_data):
    y_o = y_all[i]
    y_f = y_all[i] + alpha1*np.cos(np.angle(U_data[1,i]))
    z_o = z_all[i]
    z_f = z_all[i] + alpha2*np.sin(np.angle(U_data[1,i]))
    dat, = plt.plot([y_o,y_f],[z_o,z_f],'g')
    plt.plot(y_o,z_o,'xg')
    plt.plot(y_f,z_f,'og')

# plot uncorrelated phases for grid
for i in range(n_data,n_all):
    y_o = y_all[i]
    y_f = y_all[i] + alpha1*np.cos(np.angle(Xu_all[i]))
    z_o = z_all[i]
    z_f = z_all[i] + alpha2*np.sin(np.angle(Xu_all[i]))
    unc, = plt.plot([y_o,y_f],[z_o,z_f],'r')
    plt.plot(y_o,z_o,'xr')
    plt.plot(y_f,z_f,'or')

# plot correlated phases for all
for i in range(n_all):
    y_o = y_all[i]
    y_f = y_all[i] + alpha1*np.cos(np.angle(U_all[1,i]))
    z_o = z_all[i]
    z_f = z_all[i] + alpha2*np.sin(np.angle(U_all[1,i]))
    if i < n_data: col = 'b'
    else:          col = 'k'
    corr, = plt.plot([y_o,y_f],[z_o,z_f],col)
    plt.plot(y_o,z_o,'x'+col)
    plt.plot(y_f,z_f,'o'+col)

## plot data phases
#for i in range(n_data):
#    y_o = y_all[i]
#    y_f = y_all[i] + alpha*np.abs(U_data[1,i])*np.cos(np.angle(U_data[1,i]))
#    z_o = z_all[i]
#    z_f = z_all[i] + alpha*np.abs(U_data[1,i])*np.sin(np.angle(U_data[1,i]))
#    dat, = plt.plot([y_o,y_f],[z_o,z_f],'g')
#    plt.plot(y_o,z_o,'xg')
#    plt.plot(y_f,z_f,'og')
#
## plot uncorrelated phases for grid
#for i in range(n_data,n_all):
#    y_o = y_all[i]
#    y_f = y_all[i] + alpha*np.abs(Xu_all[i])*np.cos(np.angle(Xu_all[i]))
#    z_o = z_all[i]
#    z_f = z_all[i] + alpha*np.abs(Xu_all[i])*np.sin(np.angle(Xu_all[i]))
#    unc, = plt.plot([y_o,y_f],[z_o,z_f],'r')
#    plt.plot(y_o,z_o,'xr')
#    plt.plot(y_f,z_f,'or')
#
## plot correlated phases for all
#for i in range(n_all):
#    y_o = y_all[i]
#    y_f = y_all[i] + alpha*np.abs(U_all[1,i])*np.cos(np.angle(U_all[1,i]))
#    z_o = z_all[i]
#    z_f = z_all[i] + alpha*np.abs(U_all[1,i])*np.sin(np.angle(U_all[1,i]))
#    if i < n_data: col = 'b'
#    else:          col = 'k'
#    corr, = plt.plot([y_o,y_f],[z_o,z_f],col)
#    plt.plot(y_o,z_o,'x'+col)
#    plt.plot(y_f,z_f,'o'+col)

plt.xlim([y_grid[0]-10,y_grid[-1]+10])
plt.ylim([z_grid[0]-10,z_grid[-1]+10])

plt.figure(2)
plt.clf()
plt.plot(u_data[:,0])
plt.plot(u_all[:,0])


# initialize figure
fig = plt.figure(3)
plt.clf()

# use grid to make tick labels
sk_ang  = np.pi/8                           # plot skew angle
n_p     = 6
y_ticks = jr.grid2ticklist(y_grid)
z_ticks = jr.grid2ticklist(z_grid)
i_ts   = np.arange(n_p) * n_t/n_p           # time indices of snapshots
v_min, v_max = 4, 19
#lvls = np.linspace(v_min,v_max,15)
lvls = np.arange(v_min,v_max+1)

# loop through snapshots
for i in range(n_p):
    
    plt.subplot(1,n_p,i+1)
    plt.contourf(Y_grid,Z_grid,u_grid[:,:,i_ts[i]], \
        cmap='afmhot_r',origin='lower')
    
#    # define subplot identifier (e.g., 121)
#    i_p = int(100 + 10*n_p + i + 1)   # plot index
#    
#    # intialize sheared axes
#    ax_flt, aux_shr = jr.sheared_axes(fig, i_p, y_ticks, z_ticks,\
#                                        skew_ang = sk_ang)
#
#    # plot contour of turbulent field
#    cnt = aux_shr.contourf(Y_grid,Z_grid,u_grid[:,:,i], \
#        cmap='afmhot_r',origin='lower',levels=lvls)
#    
#    # plot grid only on first image        
#    if (i==0):
#        aux_shr.scatter(Y_grid[:],Z_grid[:],s=1,c='k',zorder=lvls.size+1)
#        
#    # add time on top of plot
#    aux_shr.text(0,2*z_grid[-1]-z_grid[-2],'t = {:.0f}'.format(i_ts[i]*dt),\
#                ha='right', fontsize='large',
#                rotation=sk_ang*180./np.pi)
    
## add colorbar to right of rightmost plot
#cbaxes = fig.add_axes([0.95, 0.1, 0.02, 0.8]) 
#cb = plt.colorbar(cnt, cax = cbaxes) 

# tighten up plot
#plt.subplots_adjust(left=0.03, right=0.92, top=1.00, bottom=0.05)

plt.figure(4)
plt.clf()
plt.plot(u_grid[2,4,:])
plt.plot(u_data[:,2])

plt.figure(5)
plt.clf()
plt.subplot(141)
plt.contourf(Y_grid,Z_grid,np.mean(u_grid,axis=2), \
        cmap='afmhot_r',origin='lower')
plt.colorbar()
plt.subplot(142)
plt.plot(np.mean(u_data,axis=0),z_data,'.')
plt.plot(p_out[0] * np.power(z_grid,p_out[1]),z_grid)
plt.subplot(143)
plt.contourf(Y_grid,Z_grid,np.std(u_grid,axis=2), \
        cmap='afmhot_r',origin='lower')        
plt.colorbar()
plt.subplot(133)
plt.plot(np.std(u_data,axis=0),z_data)
