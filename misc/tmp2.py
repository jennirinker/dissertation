"""
Compare FRF and ODE simulations for SDOF
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

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
    
    if len(u.shape) > 1:
        u = np.squeeze(u)

    x0 = [ICs(u[0],As,wns,zetas), 0.0]
    ForcedMassSpring = lambda state, t: \
        MassSpring(state,t,t_u,u,As,wns,zetas)
    
    x = odeint(ForcedMassSpring, x0, t)
    y = x[:,0]
    
    return y



# generate forcing
U,sig,tau,rho,mu   = 20.,1.,34.,0.0,np.pi
n_s = 1
n_t = 12000
dt  = 0.05
t, u = jr.generateKaimal1D(n_t,n_s,dt,U,sig,tau,rho,mu)
#t_u = np.arange(n_t)*dt
#u   = np.zeros(t_u.shape)

# switched system parameters
As = [2.,1.]
zetas = [0.01,0.01]
wns = [2*np.pi*1,2*np.pi*2]

y_sw = switchedSDOFresponse(t,u,As,zetas,wns)




# define response function
def SDOFresp(t,x):
    import numpy as np
    
    dt = t[5]/5.
    T  = t.size * dt
    df = 1./T
    
    wn   = 2*np.pi*1
    zeta = 0.01
    A    = 1.
    
    X = np.fft.rfft(x,axis=0)
    f = np.arange(X.shape[0])*df
    w = 2 * np.pi * f
    H = A/(wn**2 - np.power(w,2) + 2*zeta*wn*1j*w)
    
    Y = X * (H.reshape(X.shape[0],1)*np.ones((1,X.shape[1])))
    y = np.fft.irfft(Y,axis=0)
    
    return y
   
y_SS = SDOFresp(t,u)   
   
plt.figure(1)
plt.clf()

plt.plot(t, y_sw)
plt.plot(t,y_SS)
plt.ylim([-5,5])
plt.xlabel('TIME (sec)')
plt.ylabel('STATES')
plt.title('Mass-Spring System')
plt.legend(('$x$ (m)', '$\dot{x}$ (m/sec)'))    
    

plt.figure(2)
plt.clf()
plt.plot(t,state[:,0]-np.squeeze(y))

