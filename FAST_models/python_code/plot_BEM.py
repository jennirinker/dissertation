"""
Load an AeroDyn input file and plot wireframe model of
blade that is used in BEM calculations
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as scio
import os
from mpl_toolkits.mplot3d import Axes3D

def ReadAeroDyn(ADfname):
    """ Load aerodynamic data from AeroDyn file to visualize the BEM model
        using a wireframe model
    """
    import numpy as np
    
    # open AeroDyn file
    with open(ADfname,'r') as f:
        
        # initialize line counter, set NFoil dummy value, loop through file
        i_line = 0
        NFoil  = 100
        for line in f:
            
            # get total number of airfoils
            if (i_line == 17):
                NFoil = int(line.split()[0])
                AFNames = []
                
            # get list of airfoils w/o exension
            elif ((i_line > 17) and (i_line <= 17 + NFoil)):
                AFName = line.split()[0].split('/')[-1].split('.')[0].strip('/"')
                AFNames.append(AFName)
                
            # get total number of aerodynamic elements
            elif (i_line == 18 + NFoil):
                NADNodes = int(line.split()[0])
                ADData   = np.empty((NADNodes,5))
                i_BlD    = 0
                
            # get AD element data
            elif ((i_line > 19 + NFoil) and (i_line <= 19+NFoil+NADNodes)):
                ADData[i_BlD,:] = [float(x) for x in line.split()[:5]]
                i_BlD += 1
            
            # increment line counter
            i_line += 1
            
    fields = ['RNodes','AeroTwst','DRNodes','Chord','NFoil']
    ADDict = {}
    ADDict['ADFields'] = fields
    ADDict['Airfoils'] = AFNames
    ADDict['ADData']   = ADData
    
    return ADDict

# path to AeroDyn file
ADfname = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
            'dissertation\\FAST_models\\FAST7\\WP0750_v7\\WP0750_AD.ipt'
ADDict = ReadAeroDyn(ADfname)

# path to airfoil directory
AFdir = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation' + \
            '\\FAST_models\\python_code'

## set path to AeroDyn input file and blade geometries
#ADdir = ''
#BlGeoDir = ''
#
## read through AeroDyn file, get list of airfoil names and then a list of 
## element information
#
# initialize figure and 3D axes
fig = plt.figure(1,figsize=(8.5,3.0))
plt.clf()
axTop = plt.axes([0.12,0.60,0.80,0.35])
ax3D = Axes3D(fig,rect=[-0.05,0.02,1.05,0.46]) 

#ADData = np.array([
#            [0.5,0,1,1,0],
#            [2.0,0,2,2,0]])
#
ADData = ADDict['ADData']
Airfoils = ADDict['Airfoils']
Fields = [s.strip() for s in ADDict['ADFields']]

# loop through values in BEM schedule and incrementally add to wireframe plot
for i_el in range(len(ADData)):
    
    # extract information for that element
    RNode  = ADData[i_el,Fields.index('RNodes')]
    Twst   = ADData[i_el,Fields.index('AeroTwst')]
    DR     = ADData[i_el,Fields.index('DRNodes')]
    Chord  = ADData[i_el,Fields.index('Chord')]
    AFID = ADData[i_el,Fields.index('NFoil')]

    # get beginning and ending nodes for the element
    z_lo   = RNode - DR/2.
    z_hi   = RNode + DR/2.

    # load blade geometries
    AFName = Airfoils[int(AFID-1)]
    fpathGeo = os.path.join(AFdir,AFName+'.mat')
    BlSched  = scio.loadmat(fpathGeo,squeeze_me=True)['Geometry']
    BlGeo   = BlSched * Chord
    print(BlGeo[:4])
    poop

    # rotate to include twist
    x_Bld = BlGeo[:,0]*np.cos(Twst*np.pi/180.) - \
                BlGeo[:,1]*np.sin(Twst*np.pi/180.)
    y_Bld = BlGeo[:,0]*np.sin(Twst*np.pi/180.) - \
                BlGeo[:,1]*np.cos(Twst*np.pi/180.)
    
    # convert to a format valid for plot_surface
    X = np.vstack((x_Bld.reshape(1,x_Bld.size),
                   x_Bld.reshape(1,x_Bld.size)))
    Y = np.vstack((y_Bld.reshape(1,y_Bld.size),
                   y_Bld.reshape(1,y_Bld.size)))
    Z = np.vstack((z_lo*np.ones((1,y_Bld.size)),
                   z_hi*np.ones((1,y_Bld.size))))
                   
    # add to wireframe plot (rearrange for graphing)
    ax3D.plot_wireframe(Z,Y,X)
#    ax3D.plot_surface(X.T,Y.T,Z.T)
    
    # add to top-view plot

    axTop.plot([z_lo,z_hi,z_hi,z_lo,z_lo],[0,0,Chord,Chord,0],'b')
        
# prettify 2D plot
axTop.set_ylim([0,5]) 
axTop.set_xlim([0,np.ceil(RNode+DR)])  
axTop.spines['top'].set_visible(False)
axTop.spines['right'].set_visible(False)
axTop.get_xaxis().tick_bottom()

axTop.get_yaxis().tick_left()
           
# set axis labels
ax3D.set_xlabel('Z')
ax3D.set_ylabel('Y')
ax3D.set_zlabel('X')

#ax3D.view_init(elev=0., azim=0) 
#ax3D.w_zaxis.line.set_lw(0.)
#ax3D.set_zticks([])
#ax3D.grid('off')     
#ax3D.w_xaxis.set_pane_color((0., 0., 0., .0))
#ax3D.w_yaxis.set_pane_color((0., 0., 0., .0))
#ax3D.w_zaxis.set_pane_color((0., 0., 0., .0))
    