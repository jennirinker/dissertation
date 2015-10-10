"""
Plot lift and drag coefficients of specified airfoil
"""
import numpy as np
import scipy.io as scio
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate
#
#fname = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#            'dissertation\\FAST_models\\FAST7\\WP0750_v7\\WP0750_AD.ipt'

#with open(fname,'r') as f:
#    
#    i_line = 0
#    for line in f:
#        
#        if (i_line == 17):
#            NFoil = int(line.split()[0])
#            AFNames = []
#        elif ((i_line > 17) and (i_line <= 17 + NFoil)):
#            AFName = line.split()[0].split('/')[-1].strip('/"')
#            AFNames.append(AFName)
#        elif (i_line == 18 + NFoil):
#            NBlNodes = int(line.split()[0])
#            BlData   = np.empty((NBlNodes,5))
#            i_BlD    = 0
#        elif ((i_line > 19 + NFoil) and (i_line <= 19+NFoil+NBlNodes)):
#            BlData[i_BlD,:] = [float(x) for x in line.split()[:5]]
#            i_BlD += 1
#        
#        i_line += 1

AFs = scio.loadmat('WindPACT_geometries.mat',squeeze_me=True)
DR = 1.6
AFName = 'S825_2102'
AFGeo  = AFs[AFName]

x = AFGeo[:,0]
y = AFGeo[:,1]

X = np.vstack((x.reshape(1,x.size),x.reshape(1,x.size)))
Y = np.vstack((y.reshape(1,y.size),y.reshape(1,y.size)))
Z = np.vstack((np.zeros((1,y.size)),np.ones((1,y.size))))

#
fig = plt.figure(1)
plt.clf()
#ax = fig.gca(projection='3d')
#surf = ax.plot_surface(X, Y, Z)

#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib
#import numpy as np
#from scipy.interpolate import interp1d
#from matplotlib import cm
#from matplotlib import pyplot as plt
#step = 0.04
#maxval = 1.0
fig = plt.figure(1)
fig.clf()
ax = Axes3D(fig)  
#
#
#
#u=np.array([0,1,2,1,0,2,4,6,4,2,1])
#r=np.array([0,1,2,3,4,5,6,7,8,9,10])
#f=interp1d(r,u)
#
## walk along the circle
#p = np.linspace(0,2*np.pi,50)
#R,P = np.meshgrid(r,p)
## transform them to cartesian system
#X,Y = R*np.cos(P),R*np.sin(P)
#
#Z=f(R)
#
ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1)
#ax.set_xticks([])
#plt.show()
ax.set_ylim([-0.5,0.5])
ax.set_xlim([0,1])
