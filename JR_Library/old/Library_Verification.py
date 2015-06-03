""" Script to automatically verify subfunctions in Useful.py
"""

import numpy as np
import Library_JR as jr

pi = np.pi

# ================================ wrap =======================================
print '\nTesting \'wrap\':'
t1 = pi/2
t2 = 3*pi/2
[t1_out,t2_out]  = jr.wrap(np.array([t1,t2]))
print '  Theory:  {:.2f}  {:.2f}'.format(t1,t2-2*pi)
print '  Output:  {:.2f}  {:.2f}'.format(t1_out,t2_out)

# ========================== samplePhaseCoherence =============================
print '\nTesting \'samplePhaseCoherence\':'
t1 = 2*pi*np.random.rand(20000,1)
t2 = np.ones((20000,1))
[rho1,mu1] = jr.samplePhaseCoherence(t1)
[rho2,mu2] = jr.samplePhaseCoherence(t2)
print '  Theory: {:.2f}  {:s}  {:.2f}  {:.2f}'.format(0,'Undf',1,1)
print '  Output: {:.2f}  {:.2f}  {:.2f}  {:.2f}'.format(rho1,mu1,rho2,mu2)

# ========================== signalPhaseCoherence =============================
print '\nTesting \'signalPhaseCoherence\':'
rho_theo = 0.6
mu_theo = -2*pi/3
[t,x] = jr.generateKaimal1D(12000,1,0.05,10.,1.2,34.,rho_theo,mu_theo)
[rho_out,mu_out] = jr.signalPhaseCoherence(x)
print '  Theory:  {:.2f}  {:.2f}'.format(rho_theo,mu_theo)
print '  Output:  {:.2f}  {:.2f}'.format(rho_out,mu_out)

# ============================= uniqueComponents ==============================
print '\nTesting \'uniqueComponents\':'
n_t = 12000
n_f_out = jr.uniqueComponents(n_t)
print '  Theory:  {}'.format(6001)
print '  Output:  {}'.format(n_f_out)
    
# =================================== X2Sk ====================================
print '\nTesting \'X2Sk\':'
X = np.array([0.,1.,2.,1.])
Sk_out = jr.X2Sk(X)
print '  Theory:  {}  {}  {}'.format(0.,2.,8.)
print '  Output:  {}  {}  {}'.format(Sk_out[0],Sk_out[1],Sk_out[2])

# ================================= Sk2Xuniq ==================================
print '\nTesting \'Sk2Xuniq\':'
X_out = jr.Sk2Xuniq(Sk_out)
print '  Theory:  {}  {}  {}'.format(X[0],X[1],X[2])
print '  Output:  {}  {}  {}'.format(X_out[0],X_out[1],X_out[2])

# ================================= Xuniq2X ===================================
print '\nTesting \'Xuniq2X\':'
Xuniq = np.array([0.,1+1j,2-2j])
X_out1 = jr.Xuniq2X(Xuniq,4)
X_out2 = jr.Xuniq2X(Xuniq,5)
print '  Theory (even):  {}  {}  {}  {}'.format(Xuniq[0],Xuniq[1], \
    Xuniq[2],np.conj(Xuniq[1]))
print '  Output (even):  {}  {}  {}  {}'.format(X_out1[0],X_out1[1], \
    X_out1[2],X_out1[3])
print ''
print '  Theory (odd):   {}  {}  {}  {}  {}'.format(Xuniq[0],Xuniq[1], \
    Xuniq[2],np.conj(Xuniq[2]),np.conj(Xuniq[1]))
print '  Output (odd):   {}  {}  {}  {}  {}'.format(X_out2[0],X_out2[1], \
    X_out2[2],X_out2[3],X_out2[4])

# ================================ signal2Sk ==================================
print '\nTesting \'signal2Sk\':'
Sk_theo = np.array([0.,2.,8.])
Xmag = jr.Sk2Xuniq(Sk_theo)
X = jr.Xuniq2X(Xmag,4)
x = np.fft.ifft(X)*4
Sk_out = jr.signal2Sk(x)
print '  Theory:  {}  {}  {}'.format(Sk_theo[0],Sk_theo[1],Sk_theo[2])
print '  Output:  {}  {}  {}'.format(Sk_out[0],Sk_out[1],Sk_out[2])

# ============================== spectralScale ================================
print '\nTesting \'spectralScale\':'
Sk = np.array([0.,2.,8.])
n_t = 4
sig_theo = 2.
Xmag = jr.Sk2Xuniq(Sk_theo)
X = jr.Xuniq2X(Xmag,n_t)
x = np.real(np.fft.ifft(X))*n_t
alpha = jr.spectralScale(Sk,sig_theo,n_t)
y = x*alpha
sig_out = np.std(y,ddof=1)
print '  Theory:  {}'.format(sig_theo)
print '  Output:  {}'.format(sig_out)

# ============================= KaimalSpectrum ================================
print '\nTesting \'KaimalSpectrum\':'
f = np.array([0.2,0.3,5.])
tau = 34.
sig = 2.
S_theo = (sig**2)*(4*tau)/np.power(1.+6.*f*tau,5./3.)
S_out  = jr.KaimalSpectrum(f,tau,sig)
print '  Theory:  {:.3f}  {:.3f}  {:.3f}'.format(S_theo[0],S_theo[1],S_theo[2])
print '  Output:  {:.3f}  {:.3f}  {:.3f}'.format(S_out[0],S_out[1],S_out[2])

# ============================ calculateKaimal ================================
print '\nTesting \'calculateKaimal\':'
n_t = 12000
n_m = 1
dt = 0.05
U_theo = 10.
sig_theo = 1.4
tau_theo = 31.
rho_theo = 0.2
mu_theo = pi/2.
[t,x]=jr.generateKaimal1D(n_t,n_m,dt,U_theo,sig_theo,tau_theo,rho_theo,mu_theo)
tau_out = jr.calculateKaimal(x,dt)
print '  Theory:  {:.2f}'.format(tau_theo)
print '  Output:  {:.2f}'.format(tau_out)

# ============================ generateKaimal1D ===============================
print '\nTesting \'generateKaimal1D\':'
U_out = np.mean(x)
sig_out = np.std(x,ddof=1)
[rho_out,mu_out] = jr.signalPhaseCoherence(x)
print '  Theory:  {:.2f}  {:.2f}'.format(U_theo,sig_theo)
print '  Output:  {:.2f}  {:.2f}'.format(U_out,sig_out)
print '  Theory:  {:.2f}'.format(tau_theo)
print '  Output:  {:.2f}'.format(tau_out)
print '  Theory:  {:.2f}  {:.2f}'.format(rho_theo,mu_theo)
print '  Output:  {:.2f}  {:.2f}'.format(rho_out,mu_out)


# ========================== wrappedCauchySample ==============================
print '\nTesting \'wrappedCauchySample\':'
n_t = 12000
n_m = 1
mu_theo  = pi/6
rho_theo = 0.4
theta = jr.wrappedCauchySample(n_t,n_m,rho_theo,mu_theo)
[rho_out,mu_out] = jr.samplePhaseCoherence(theta)
print '  Theory:  {:.2f}  {:.2f}'.format(rho_theo,mu_theo)
print '  Output:  {:.2f}  {:.2f}'.format(rho_out,mu_out)

print '\nVerification complete'




