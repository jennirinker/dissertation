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

# ============= inputs =============

# define locations of data points
y_data = np.array([0.,0.,0.,0.,0.,0.])
z_data = np.array([15.,30,50,76,100,131.])
n_data  = y_data.size

# define desired grid
y_grid = np.linspace(-15,15,9)
z_grid = np.linspace(20,160,10)

# ============= create artificial time series =============
n_t, dt = 12000, 0.05
zhub, Vhub, turbc = 90., 8., 'A'
U_z  = jr.IEC_VelProfile(z_data,zhub,Vhub)
Iref = jr.IEC_Iref(turbc)
sig1 = jr.IEC_Sigma1(Iref,Vhub)
L = 8.1 * 42
u_data = np.empty((n_t,n_data))
for i in range(n_data):
    u_data[:,i] = np.squeeze(jr.generateKaimal1D(n_t,1,dt, \
        U_z[i],sig1,L/U_z[i],0.0,0.0)[1])
        
# ============= grid, delta_r, magnitudes, velocity profile =============
        
# useful numbers for matrix indexing
n_y    = y_grid.size
n_z    = z_grid.size
n_grid  = n_y*n_z
n_all   = n_data + n_grid
n_f     = jr.uniqueComponents(n_t)
df      = 1./(n_t*dt)
fs      = np.arange(n_f)*df

# get arrays of all points, data and grid
Y_grid, Z_grid = np.meshgrid(y_grid,z_grid)
y_all = np.concatenate((y_data,Y_grid.reshape(n_grid)))
z_all = np.concatenate((z_data,Z_grid.reshape(n_grid)))

# calculate matrix of distances
DR = np.empty((n_all,n_all))
for i in range(n_all):
    for j in range(i,n_all):
        dy  = y_all[i] - y_all[j]
        dz  = z_all[i] - z_all[j]
        DR[i,j] = np.sqrt(dy**2 + dz**2)
        DR[j,i] = np.sqrt(dy**2 + dz**2)
        
# calculate magnitudes from data, and assign magnitudes to grid FFT array
U_data = np.fft.rfft(u_data,axis=0)/n_t
Umags_all = np.empty((n_f,n_all))
Umags_all[:,:n_data] = np.abs(U_data)
for i in range(n_data,n_all):
    
    z_p = z_all[i]                          # grid point height
    dzs = np.abs(z_data-z_p)                # diffs between msmts and grid pts
    i_closest = np.where(dzs == dzs.min())  # closest measurement
    S_closest = np.power( \
                Umags_all[:,i_closest],2)   # S(f) closest point(s)
    U_p       = np.squeeze(np.sqrt( \
                np.mean(S_closest,axis=1))) # mean-power magnitudes
    Umags_all[:,i] = U_p                    # save magnitude
    
# fit velocity profile to data, save results in Umag_all
U_z = np.mean(u_data,axis=0)
fitfunc = lambda p, x: p[0] * x ** (p[1])
errfunc = lambda p, x, y: (y - fitfunc(p, x))
p_out = optimize.leastsq(errfunc, [10.,0.2],args=(z_data, U_z))[0]
Umags_all[0,n_data:] = p_out[0] * np.power(z_all[n_data:],p_out[1])
    
# ============= spatially correlate phases, create field =============
U_all = np.empty((n_f,n_all),dtype='complex')
U_all[:,:n_data] = U_data
U_all[0,n_data:] = Umags_all[0,n_data:]
for i in range(1,n_f):
    
    # get Cholesky decomposition
    f = fs[i]
    Coh = jr.IEC_SpatialCoherence(zhub,Vhub,DR,f)
    C = np.linalg.cholesky(Coh)
    
    # solve for un-spatially correlated data Fourier coefficients
    a = C[:n_data,:n_data]
    b = U_data[i,:].reshape(n_data,1)
    Xu_data = np.linalg.solve(a,b)
    
    # Fourier component for grid
    phis_grid  = np.random.rand(n_grid,1) * 2 * np.pi
    Umags_grid = Umags_all[i,n_data:].reshape(n_grid,1)
    Xu_grid    =  Umags_grid * np.exp(1j*phis_grid)
    
    # create entire vector
    Xu_all  = np.concatenate((Xu_data,Xu_grid),axis=0)
    
    # correlate
    Xs_all = np.dot(C,Xu_all)
    
    # save in array
    U_all[i,:] = np.squeeze(Xs_all)

# take IFT to get time series
u_all = np.fft.irfft(U_all,axis=0)*n_t
u_grid = np.empty((n_z,n_y,n_t))
i_tot = n_data
for i in range(n_y):
    for j in range(n_z):
        u_grid[j,i,:] = u_all[:,i_tot]
        i_tot += 1

# %%  ============= test plot =============
#
## create plot
#plt.figure(1,figsize=(6,6))
#plt.clf()
#ax = plt.axes([0.12,0.12,0.81,0.66])
#alpha = 5
#
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
#
#plt.xlim([y_grid[0]-10,y_grid[-1]+10])
#plt.ylim([z_grid[0]-10,z_grid[-1]+10])

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

### *************** NOOOO GOOD! AXES LOOK LIGHT THEY MIGHT BE REVERESED