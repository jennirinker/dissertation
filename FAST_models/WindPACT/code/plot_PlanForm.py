"""
Load an AeroDyn input file and plot wireframe model of
blade that is used in BEM calculations
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as scio
import os
from mpl_toolkits.mplot3d import Axes3D

## set path to AeroDyn input file and blade geometries
#ADdir = ''
#BlGeoDir = ''
#
## read through AeroDyn file, get list of airfoil names and then a list of 
## element information
#

NSplines = 5

# initialize figure and 3D axes
fig = plt.figure(1,figsize=(8.5,6.0))
plt.clf()
ax3D = Axes3D(fig,rect=[0.02,0.02,0.96,0.46]) 

BlData = np.array([
            [0.05,0.5,1],
            [0.07,0.3,1],
            [0.25,3.0,2],
            [0.75,1.0,2],
            [1.00,0.75,2]])

BlSched1 = np.array([[1,0],[0,1],[-1,0],[0,-1],[1,0]])
BlSched2 = np.array([[0.75,0],[0,0.4],[-0.25,0],[0,-0.4],[0.75,0]])

X_wire = np.empty((len(BlData),NSplines))
Y_wire = np.empty((len(BlData),NSplines))
Z_wire = np.empty((len(BlData),NSplines))

# loop through values in BEM schedule and incrementally add to wireframe plot
for i_el in range(len(BlData)):
    
    # extract information for that element
    z     = BlData[i_el,0]
    Chord = BlData[i_el,1]
    AFID  = BlData[i_el,2]
    
    # load blade geometries
    if AFID == 1:
        BlSched = BlSched1
    else:
        BlSched = BlSched2
    BlGeo   = BlSched * Chord
    
    # convert to a format valid for plot_surface
    x = BlGeo[:,0]
    y = BlGeo[:,1]
    
    # resample to an even number of points
    x_int = x
    y_int = y
    
    # add to wireframe array
    X_wire[i_el,:] = x_int
    Y_wire[i_el,:] = y_int
    Z_wire[i_el,:] = z
                   
    # add to wireframe plot (rearrange for graphing)
ax3D.plot_wireframe(Z_wire,X_wire,Y_wire)
           
# set axis labels
ax3D.set_xlabel('Z')
ax3D.set_ylabel('X')
#ax3D.set_zlabel('Y')

ax3D.view_init(elev=60., azim=-70) 
ax3D.w_zaxis.line.set_lw(0.)
ax3D.set_zticks([])
ax3D.grid('off')     
ax3D.w_xaxis.set_pane_color((0., 0., 0., .0))
ax3D.w_yaxis.set_pane_color((0., 0., 0., .0))
ax3D.w_zaxis.set_pane_color((0., 0., 0., .0))
    