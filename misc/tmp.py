
"""
Testing expansion of time series points to new grid
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np

# ============= inputs =============

# define locations of data points
#y_data = np.array([0.,0.,0.,0.,0.,0.])
#z_data = np.array([15.,30,50,76,100,131.])
y_data = np.array([0.,0.])
z_data = np.array([80.,90])
n_data  = y_data.size

# define desired grid
y_grid = np.linspace(-15,15,3)
z_grid = np.linspace(75,105,4)

# ============= create artificial time series =============
n_t, dt = 12000, 0.05
zhub, Vhub, turbc = 90., 8., 'A'
U_z  = jr.IEC_VelProfile(z_data,zhub,Vhub)
Iref = jr.IEC_Iref(turbc)
sig1 = jr.IEC_Sigma1(Iref,Vhub)
L = 8.1 * 42
u_all = np.empty((n_t,n_data))
for i in range(n_data):
    u_all[:,i] = np.squeeze(jr.generateKaimal1D(n_t,1,dt, \
        U_z[i],sig1,L/U_z[i],0.0,0.0)[1])

# ============= create turbulent field =============

## get Fourier components of turbulent data
#U_all = np.fft.rfft(u_all,axis=0)/n_t
#
#for i in range(n_grid):
    #ENDED HERE!!!

## useful numbers for matrix indexing
n_y    = y_grid.size
n_z    = z_grid.size
n_grid  = n_y*n_z
n_all   = n_data + n_grid

# get arrays of all points, data and grid
Z_grid, Y_grid = np.meshgrid(z_grid,y_grid)
y_all = np.concatenate((y_data,Y_grid.reshape(n_grid)))
z_all = np.concatenate((z_data,Z_grid.reshape(n_grid)))

## calculate matrix of distances
#DR = np.empty((n_all,n_all))
#for i in range(n_all):
#    for j in range(i,n_all):
#        dy  = y_all[i] - y_all[j]
#        dz  = z_all[i] - z_all[j]
#        DR[i,j] = np.sqrt(dy**2 + dz**2)
#        DR[j,i] = np.sqrt(dy**2 + dz**2)
#
## create field for each frequency
#n_f = jr.uniqueComponents(n_t)
#df  = 1./(n_t*dt)
##fs  = np.arange(1,n_f)*df
#fs = [0.,0.1]
#for i in range(1,2):
#    
#    # get Cholesky decomposition
#    f = fs[i]
#    Coh = jr.IEC_SpatialCoherence(zhub,Vhub,DR,f)
#    C = np.linalg.cholesky(Coh)
#    
#    # solve for un-spatially correlated data Fourier coefficients
#    a = C[:n_data,:n_data]
#    b = U_all[i,:].reshape(n_data,1)
#    Xu_data = np.linalg.solve(a,b)
#    
#    # make up magnitudes for other Fourier coefficients
#    Xu_grid = 1.1*np.exp(1j*2*np.pi*np.random.rand(n_grid,1))
#    
#    # create entire vector
#    Xu_all = np.concatenate((Xu_data,Xu_grid),axis=0)
#    
#    # correlate
#    Xs_all = np.dot(C,Xu_all)
#    

# construct matrix of spatial coherences
Coh = np.empty((n_all,n_all))
for i in range(n_all):
    for j in range(i,n_all):
        dy  = y_all[i] - y_all[j]
        dz  = z_all[i] - z_all[j]
        dr  = np.sqrt(dy**2 + dz**2)
        coh = jr.IEC_SpatialCoherence(zhub,Vhub,dr,f)
        Coh[i,j] = coh
        Coh[j,i] = coh

# get Cholesky matrix
C = np.linalg.cholesky(Coh)

# make up known Fourier coefficients
Xs_data = np.array([[1.1*np.exp(1j*np.pi/4)],[1.15*np.exp(1j*5*np.pi/16)]])
#Xs_data = np.array([[1.1*np.exp(1j*np.pi/4)]])

# solve for USC data Fourier coefficients
a = C[:n_data,:n_data]
b = Xs_data
Xu_data = np.linalg.solve(a,b)

# make up magnitudes for other Fourier coefficients
Xu_grid = 1.1*np.exp(1j*2*np.pi*np.random.rand(n_grid,1))

# create entire vector
Xu_all = np.concatenate((Xu_data,Xu_grid),axis=0)

# correlate
Xs_all = np.dot(C,Xu_all)

# create plot
plt.figure(1,figsize=(6,6))
plt.clf()
ax = plt.axes([0.12,0.12,0.81,0.66])
alpha = 2.

# plot data phases
for i in range(n_data):
    y_o = y_all[i]
    y_f = y_all[i] + alpha*np.abs(Xs_data[i])*np.cos(np.angle(Xs_data[i]))
    z_o = z_all[i]
    z_f = z_all[i] + alpha*np.abs(Xs_data[i])*np.sin(np.angle(Xs_data[i]))
    dat, = plt.plot([y_o,y_f],[z_o,z_f],'b')
    plt.plot(y_o,z_o,'xb')
    plt.plot(y_f,z_f,'ob')

# plot uncorrelated phases for grid
for i in range(n_data,n_all):
    y_o = y_all[i]
    y_f = y_all[i] + alpha*np.abs(Xu_all[i])*np.cos(np.angle(Xu_all[i]))
    z_o = z_all[i]
    z_f = z_all[i] + alpha*np.abs(Xu_all[i])*np.sin(np.angle(Xu_all[i]))
    unc, = plt.plot([y_o,y_f],[z_o,z_f],'r')
    plt.plot(y_o,z_o,'xr')
    plt.plot(y_f,z_f,'or')

# plot correlated phases for all
for i in range(n_all):
    y_o = y_all[i]
    y_f = y_all[i] + alpha*np.abs(Xs_all[i])*np.cos(np.angle(Xs_all[i]))
    z_o = z_all[i]
    z_f = z_all[i] + alpha*np.abs(Xs_all[i])*np.sin(np.angle(Xs_all[i]))
    if i < n_data: col = 'b'
    else:          col = 'k'
    corr, = plt.plot([y_o,y_f],[z_o,z_f],col)
    plt.plot(y_o,z_o,'x'+col)
    plt.plot(y_f,z_f,'o'+col)

plt.xlim([y_grid[0]-5,y_grid[-1]+5])
plt.ylim([z_grid[0]-5,z_grid[-1]+5])
plt.xlabel('Lateral (y)')
plt.ylabel('Vertical (z)')

plt.legend([dat,unc,corr],['Data','Pre-Spatial Correlation', \
    'Post Spatial Correlation'], 
    bbox_to_anchor=(1.0, 1.03), loc=4, borderaxespad=0.)

