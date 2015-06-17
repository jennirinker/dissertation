"""
A demonstration of the effect of spatial coherence on temporal coherence for
points with the same PDDs.

Jenni Rinker, 18-May-2015.
"""
import numpy as np
import sys
sys.path.append('..')
import JR_Library.main as jr
import matplotlib.pyplot as plt

# Round 1 calculations
n_t = 12000.
dt  = 0.05
U   = 10.
sig = 1.5
tau = 34.
rho = 0.3
mu  = np.pi

rsep = 1.0

N_mc = 500

rho0s_unc = np.empty(N_mc)
rho1s_unc = np.empty(N_mc)
rho1s_cor = np.empty(N_mc)

# run monte carlo
for i_mc in range(N_mc):

    # =========== create two un-spatially correlated simulations ===========
    t, x = jr.generateKaimal1D(n_t,2,dt,U,sig,tau,rho,mu)
    
    # ========= spatially correlate second simulation with the first =========
    n_f = jr.uniqueComponents(n_t)
    f   = np.arange(n_f)*(1/(n_t*dt))
    
    # calculate coherences at different frequencies
    Coh = jr.SpatialCoherence(90.,U,rsep,f)
    
    # get phasors of uncorrelated simulations
    X0_unc = np.fft.fft(x[:,0])[:n_f]/n_t
    X0p_unc = X0_unc/np.abs(X0_unc)
    X1_unc = np.fft.fft(x[:,1])[:n_f]/n_t
    X1m_unc = np.abs(X1_unc)
    X1p_unc = X1_unc/X1m_unc
    
    # loop through frequencies to correlate second record
    X1p_cor = np.empty(X1p_unc.shape,dtype=complex)
    X1p_cor[0] = X1p_unc[0]
    for i in range(1,n_f):
        X1p_cor[i] = Coh[i]*X0p_unc[i] + np.sqrt(1-Coh[i]**2)*X1p_unc[i]
        
    # get time history
    X1_cor  = jr.Xuniq2X(X1m_unc*X1p_cor,n_t)
    x1_cor = np.fft.ifft(X1_cor)*n_t
    
    # calculate rho, mu of newly correlated record
    rho0_unc, mu0_unc = jr.signalPhaseCoherence(x[:,0])
    rho1_unc, mu1_unc = jr.signalPhaseCoherence(x[:,1])
    rho1_cor, mu1_cor = jr.signalPhaseCoherence(x1_cor)
    
    # save concentration parameters
    rho0s_unc[i_mc] = rho0_unc
    rho1s_unc[i_mc] = rho1_unc
    rho1s_cor[i_mc] = rho1_cor
    
    # print update
    if (not i_mc % 25):
        print(i_mc)
    

#print('{:.4f}  {:.2f}'.format(rho0_unc,mu0_unc))
#print('{:.4f}  {:.2f}'.format(rho1_unc,mu1_unc))
#print('{:.4f}  {:.2f}'.format(rho1_cor,mu1_cor))

print('{:.4f}'.format(np.sum(rho0s_unc > 0.03)/float(N_mc)*100))
print('{:.4f}'.format(np.sum(rho1s_unc > 0.03)/float(N_mc)*100))
print('{:.4f}'.format(np.sum(rho1s_cor > 0.03)/float(N_mc)*100))

# comparison of time histories before and after correlation
plt.figure(1)
plt.clf()
plt.subplot(3,1,1)
plt.plot(t,x[:,0])
plt.ylabel('Point 0, Unc.')
plt.subplot(3,1,2)
plt.plot(t,x[:,1])
plt.ylabel('Point 1, Unc.')
plt.subplot(3,1,3)
plt.plot(t,x1_cor)
plt.xlabel('Time')
plt.ylabel('Point 1, Cor..')

# plot of spatial coherence function
plt.figure(2)
plt.clf()
plt.plot(f,Coh)
plt.ylim([0,1])
plt.xlabel('Frequency')
plt.ylabel('Spatial Coherence')

# histogram comparison
plt.figure(3)
plt.clf()
plt.hist(rho0s_unc,normed=True,histtype='step')
plt.hist(rho1s_unc,normed=True,histtype='step')
plt.hist(rho1s_cor,normed=True,histtype='step')
plt.legend(labels=[r'$\rho_{0,unc}$',r'$\rho_{1,unc}$',r'$\rho_{1,cor}$'])
plt.xlabel(r'$\rho$')
plt.ylabel('Normalized Histogram')

plt.show()