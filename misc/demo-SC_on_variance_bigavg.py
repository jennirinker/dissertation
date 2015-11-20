"""
investigating how spatial coherence affects variance
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
import matplotlib.pyplot as plt

# inputs
zg = np.arange(74,94.1,5)
yg = np.arange(-10,10.1,5)
#zg = np.arange(74,94.1,5)
#yg = np.array([-5,0,5])
T,dt = 630.,0.05
rho,mu = 0.2,np.pi
tau,sig = 95.2,4.2
URef = 10.5
ZRef = 100.
ZHub = 84.

nbins = 50
spb = 100
n_s = nbins*spb

# intermediate parameters
n_p = zg.size*yg.size
n_t = T/dt
n_f = jr.uniqueComponents(n_t)
df  = 1./T
fs  = np.arange(n_f)*df
UHub = URef*(ZHub/ZRef)**0.2
t = np.arange(n_t)*dt

# calculate DR
Yg,Zg = np.meshgrid(yg,zg)
yp,zp = Yg.reshape(Yg.size),Zg.reshape(Zg.size)
DR    = np.empty((n_p,n_p))
for i1 in range(n_p):
    for i2 in range(i1,n_p):
        dy = yp[i1]-yp[i2]
        dz = zp[i1]-zp[i2]
        dr = np.sqrt(dy**2+dz**2)
        DR[i1,i2] = dr
        DR[i2,i1] = dr
 

# calculate C matrix
print('\nCalculating Cholecky decompositions...')
Cs = np.empty((n_p,n_p,n_f))
for i_f in range(1,n_f):
    f = fs[i_f]
    
    # calculate Coh, C
    Coh = jr.IEC_SpatialCoherence(ZHub,UHub,DR,f)
    Cs[:,:,i_f] = np.linalg.cholesky(Coh)

        
# initialize sigma and rho
sigsX = np.empty((n_p,n_s))
rhosX = np.empty((n_p,n_s))
sigsY = np.empty((n_p,n_s))
rhosY = np.empty((n_p,n_s))        
        
for ibin in range(nbins):
    print('Processing bin {:d}...'.format(ibin))

    # create phases
#    DPhi = jr.wrappedCauchySample((n_p,spb,n_f-1),rho,mu)
#    Phi = np.zeros((n_p,spb,n_f))
#    Phi[:,:,1:] = np.cumsum(DPhi,axis=2)
    DPhi = jr.wrappedCauchySample((n_p,spb,n_f-1),rho,mu)
    DPhi[:,:,0] = np.random.rand(n_p,spb)*2*np.pi
    Phi = np.zeros((n_p,spb,n_f))
    Phi[:,:,1:] = np.cumsum(DPhi,axis=2)
    
    # get magnitudes
    S = jr.KaimalSpectrum(fs,tau,sig)
    Sk = S*df
    alpha = jr.spectralScale(Sk,sig,n_t)
    Xmag1 = jr.Sk2Xuniq(Sk)*alpha
    Xmag1[0] = UHub
    Xmag = np.ones((n_p,spb,n_f))*Xmag1.reshape(1,Xmag1.size)
    
    # define X
    X = Xmag * np.exp(1j*Phi)
    
    # initialize Y
    Y = np.empty(X.shape,dtype=complex)
    Y[:,:,0] = X[:,:,0]
    
    for i_f in range(1,n_f):
        f = fs[i_f]
        
        # extract C
        C = Cs[:,:,i_f]
        
        # set Y for i_f
        Y[:,:,i_f] = np.dot(C,X[:,:,i_f])     
    
    # calculate Xall,Yall, take IFFT of both
    xs = np.fft.irfft(X,axis=2)*n_t
    ys = np.fft.irfft(Y,axis=2)*n_t
    
    # for each point    
    for i_p in range(n_p):
    
        # for each simulation
        for i_s in range(spb):
            
            x = xs[i_p,i_s,:]
            y = ys[i_p,i_s,:]
    
            # caluclate and save sig and rho
            sigsX[i_p,i_s+ibin*spb] = np.std(x)
            rhosX[i_p,i_s+ibin*spb] = jr.signalPhaseCoherence(x)[0]
            sigsY[i_p,i_s+ibin*spb] = np.std(y)
            rhosY[i_p,i_s+ibin*spb] = jr.signalPhaseCoherence(y)[0]

#%% ================== do shit on the results ==================

# MEANS AND STANDARD DEVIATIONS
mnssigsX = np.mean(sigsX,axis=1)
mnssigsY = np.mean(sigsY,axis=1)
mnsrhosX = np.mean(rhosX,axis=1)
mnsrhosY = np.mean(rhosY,axis=1)
prcerrsigX = (mnssigsX - sig)/sig * 100
prcerrsigY = (mnssigsY - sig)/sig * 100
prcerrrhoX = (mnsrhosX - rho)/rho * 100
prcerrrhoY = (mnsrhosY - rho)/rho * 100


#print('\nMean standard deviation by point:')
#print('---------------------------------')
#print('x = {:5.3f} y = {:5.3f}'.format(mnssigsX[0],mnssigsY[0]))
#for ip in range(1,n_p):
#    print('    {:5.3f}     {:5.3f}'.format(mnssigsX[ip],mnssigsY[ip]))
#
#print('\nMean concentration parameter by point:')
#print('--------------------------------------')
#print('x = {:5.3f} y = {:5.3f}'.format(mnsrhosX[0],mnsrhosY[0]))
#for ip in range(1,n_p):
#    print('    {:5.3f}     {:5.3f}'.format(mnsrhosX[ip],mnsrhosY[ip]))
          
#sigXplot = mnssigsX
#sigYplot = mnssigsY
#rhoXplot = mnsrhosX
#rhoYplot = mnsrhosY
sigXplot = prcerrsigX
sigYplot = prcerrsigY
rhoXplot = prcerrrhoX
rhoYplot = prcerrrhoY
          
# plot means by point
plt.figure(2,figsize=(6.5,5.))
plt.clf()


#     sig - by point index
plt.subplot(3,1,1)
xplot = np.arange(n_p)+1
plt.scatter(xplot,sigYplot)
plt.grid('on')
plt.xlabel('Turb. Std. Dev. [m$^2$/s$^2$]')
plt.title('Mean of $\sigma_u$ by point index (after SC)')

#     sig - by distance to first point
plt.subplot(3,1,2)
xplot = np.sqrt((yp-yp[0])**2+(zp-zp[0])**2)
plt.scatter(xplot,sigYplot)
plt.grid('on')
plt.xlabel('Turb. Std. Dev. [m$^2$/s$^2$]')
plt.title('Mean of $\sigma_u$ by distance to Point 1 (after SC)')

#     rho
plt.subplot(3,1,3)
plt.scatter(np.arange(n_p)+1,rhoYplot)
plt.grid('on')
plt.xlabel('Concentration parameter [-]')
plt.title('Mean of rho by point index (after SC)')
              
plt.tight_layout()

# ============ X plot means by point
plt.figure(5,figsize=(6.5,5.))
plt.clf()


#     sig - by point index
plt.subplot(3,1,1)
xplot = np.arange(n_p)+1
plt.scatter(xplot,sigXplot)
plt.grid('on')
plt.xlabel('Turb. Std. Dev. [m$^2$/s$^2$]')
plt.title('Mean of $\sigma_u$ by point index (before SC)')

#     sig - by distance to first point
plt.subplot(3,1,2)
xplot = np.sqrt((yp-yp[0])**2+(zp-zp[0])**2)
plt.scatter(xplot,sigXplot)
plt.grid('on')
plt.xlabel('Turb. Std. Dev. [m$^2$/s$^2$]')
plt.title('Mean of $\sigma_u$ by distance to Point 1 (before SC)')

#     rho
plt.subplot(3,1,3)
plt.scatter(np.arange(n_p)+1,rhoXplot)
plt.grid('on')
plt.xlabel('Concentration parameter [-]')
plt.title('Mean of rho by point index (before SC)')
              
plt.tight_layout()

# pseudocolor
plt.figure(4,figsize=(4.5,4.0))
plt.clf()
SigsY = (sigYplot.reshape(Yg.shape))
plt.pcolor(Yg,Zg,SigsY,cmap='Reds')
plt.colorbar()
plt.axis('tight')
plt.xlabel('Lateral Direction [m]')
plt.ylabel('Vertical Direction [m]')
plt.title('Pseudocolor of Mean of $\sigma_u$')
plt.tight_layout()

# HISTOGRAMS
plt.figure(3,figsize=(6.5,5.))
plt.clf()
nplot = 3

#     sig
plt.subplot(2,1,1)
ip_plot = [int(i) for i in np.linspace(2,n_p-1,nplot)]
for iplot in range(nplot):
    plt.hist(sigsY[ip_plot[iplot]],
             histtype='step',label='Point {:d}'.format(ip_plot[iplot]))
plt.legend()
plt.xlabel('Turb. Std. Dev. [m$^2$/s$^2$]')
plt.title('Distribution of sig')

#     rho
plt.subplot(2,1,2)
ip_plot = [int(i) for i in np.linspace(1,n_p-1,nplot)]
for iplot in range(nplot):
    plt.hist(rhosY[ip_plot[iplot]],
             histtype='step',label='Point {:d}'.format(ip_plot[iplot]))
plt.legend()
plt.xlabel('Concentration parameter [-]')
plt.title('Distribution of rho')

plt.tight_layout()

# ================== spectra ==================
plt.figure(6,figsize=(4.5,3.5))
plt.clf()
X_avg = np.mean(X,axis=1)
Y_avg = np.mean(Y,axis=1)
plt.loglog(fs[1:],Xmag1[1:])
plt.loglog(fs[1:],np.abs(X_avg[0,1:]))
plt.loglog(fs[1:],np.abs(Y_avg[-1,1:]))

