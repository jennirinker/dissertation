"""
messing around with RSMs
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
    

def ICs(u0,As,wns,zetas):
    if u0 < 15:
        wn   = wns[0]
        zeta = zetas[0]
        A = As[0]
    else:
        wn   = wns[1]
        zeta = zetas[1]
        A = As[1]
        
    IC = u0*A/(wn**2)
    
    return IC

def MassSpring(state,t,t_u,u,As,wns,zetas):
    # unpack the state vector
    x = state[0]
    xd = state[1]
        
    # interpolated forcing
    f = np.interp(t,t_u,u)
    
    # these are our constants
    if f < 15:
        wn   = wns[0]
        zeta = zetas[0]
        A = As[0]
    else:
        wn   = wns[1]
        zeta = zetas[1]
        A = As[1]

    # compute acceleration xdd
    xdd = (-wn**2)*x + (-2*zeta*wn)*xd + A*f

    # return the two state derivatives
    return [xd, xdd]

def switchedSDOFresponse(t_u,u,As,zetas,wns):
    from scipy.integrate import odeint
    
    if len(u.shape) > 1:
        u = np.squeeze(u)

    x0 = [ICs(u[0],As,wns,zetas), 0.0]
    ForcedMassSpring = lambda state, t: \
        MassSpring(state,t,t_u,u,As,wns,zetas)
    
    x = odeint(ForcedMassSpring, x0, t_u)
    y = x[:,0]
    
    return y
    
# define inputs for system
L   = 340.
n_s = 3
n_t = 3000
dt  = 0.10
mu = np.pi
Uvec  = np.array([5,10,15,20,25])
Ivec = np.array([0,0.05,0.10,0.15,0.20])
p_i = [6,6]
rho = 0.0

# switched system parameters
As = [1.,1.5]
zetas = [0.1,0.1]
wns = [2*np.pi*1,2*np.pi*1.5]
##I = 0.15
##rhovec = np.array([0,0.1,0.3,0.5,0.7])

## ============ THIS TAKES LONG -- COMMENT OUT ONCE GENERATED
# generate responses
#n_all = Uvec.size * Ivec.size * n_s
#x = np.empty((n_all,2))
#y_resp = np.empty(n_all)
#i_all = 0
#for i_U in range(Uvec.size):        # y for plotting
#    U = Uvec[i_U]
#    for i_I in range(Ivec.size):    # x for plotting
#        I = Ivec[i_I]
##        rho = Ivec[i_I]
#        sig = I*U
#        tau = L/U
#        t, u = jr.generateKaimal1D(n_t,n_s,dt,U,sig,tau,rho,mu)
#        y = np.empty((n_t,n_s))
#        for i_s in range(n_s):
#            y[:,i_s] = switchedSDOFresponse(t,u[:,i_s],As,zetas,wns)
#        y_max = np.max(y,axis=0)
#        
#        x[i_all:i_all+n_s,0] = U
#        x[i_all:i_all+n_s,1] = I
##        x[i_all:i_all+n_s,1] = rho
#        y_resp[i_all:i_all+n_s] = y_max
#        
#        i_all += n_s
#        
#        print(i_all/float(n_all))
    
# fit polynomial surface
coeffs, ps = jr.polyregression(x,y_resp,p_i)  

# calculate surface values at data points
Iplot = np.linspace(0.05,0.20,5)
Uplot = np.linspace(Uvec[0],Uvec[-1],5)
#Iplot = np.linspace(Ivec[0],Ivec[-1],11)
#Uplot = np.linspace(Uvec[0],Uvec[-1],11)
Iplot_arr, Uplot_arr = np.meshgrid(Iplot,Uplot)
xplot = np.hstack((Uplot_arr.reshape((Uplot_arr.size,1)),
                   Iplot_arr.reshape((Iplot_arr.size,1))))
A = jr.myvander(xplot,ps)                
#A = myvander(x,ps)
y_surf = np.dot(A,coeffs)

# create results for plot
#X, Y = np.meshgrid(Ivec,Uvec)
#Z    = y_surf[0:-1:n_s].reshape((Uvec.size,Ivec.size))
X = Iplot_arr
Y = Uplot_arr
Z    = y_surf.reshape((Uplot.size,Iplot.size))
    
# plot results
fig = plt.figure(1)
plt.clf()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x[:,1],x[:,0],y_resp)

#Uplot = np.linspace(0,25,20)
#Iplot = np.linspace(0,0.20,20)

#X,Y = np.meshgrid(Uvec,rhovec)

#ax.plot_wireframe(Ugrid,Igrid,y_surf.T)
ax.plot_wireframe(X,Y,Z)
plt.xlabel('Turbulence Intensity [-]')
plt.ylabel('Mean Wind Speed [m/s]')
plt.title('SDOF Response to Wind and Meta-Model')

#fig = plt.figure(1)


