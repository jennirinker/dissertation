
"""
Testing expansion of time series points to new grid
"""
import numpy as np

y_data = np.array([0,0,0,0,0,0.])
z_data = np.array([15,30,50,76.,100,131.])

y_grid = np.linspace(-40,40,5)
z_grid = np.linspace(20,160,7)
n_y    = y_grid.size
n_z    = z_grid.size

n_data  = y_data.size
n_grid  = n_y*n_z

Y_grid, Z_grid = np.meshgrid(y_grid,z_grid)

y_all = np.concatenate((y_data,Y_grid.reshape(n_grid)))
z_all = np.concatenate((z_data,Z_grid.reshape(n_grid)))

plt.scatter(y_all,z_all)

#coh = 0.5
#C = np.array([[1.,0.],[coh,np.sqrt(1-coh**2)]])
#
#X1u = 1.1*np.exp(1j*np.pi/4)
#
#U1u = X1u
#U2u = np.abs(U1u) * np.exp(1j*np.random.rand()*2*np.pi)
#
#X2u = np.dot(C[1,:],np.array([[U1u],[U2u]]))