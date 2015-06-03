""" Subfunctions used in simulating turbulent wind in Python.

Jenni Rinker, 26-Mar-2015.
"""

def wrappedCauchySample(n_t,n_m,rho,mu):
    """ Random numbers from a Wrapped Cauchy distribution (n_t x n_m) with 
        location parameter mu and concentration parameter rho.
        
        Reference: Statistical Analysis of Circular Data, Fisher, Sec. 3.3.4
        Modified to correctly implement arccos.
        
        Args:
            n_t (int): sample size
            n_m (int): number of samples
            rho (float): concentration parameter
            mu (float): location parameter
    
        Returns:
            theta (numpy array): sample of angles
    """
    import numpy as np
    import misc
    
    if ((rho < 0) or (rho>1)):
        print 'rho must be between 0 and 1'
        return []
    
    U = np.random.rand(n_t,n_m)                     # uniform random numbers
    
    V = np.cos(2.*np.pi*U)                          # See lines before Eq. 3.28
    c = 2.*rho/(1+(rho**2))                         # See lines before Eq. 3.28
    
    B = 2*np.round(np.random.rand(n_t,n_m)) - 1     # boolean RV
                
    theta = misc.wrap(np.multiply(B, (np.arccos(np.divide(V+c, \
        1+c*V)))) + mu)                             # sample of angles
    
    return theta

def generateKaimal1D(n_t,n_m,dt,U,sig,tau,rho,mu):
    """ Single-point turbulent time history with the Kaimal spectrum and the 
        given parameters.

        Args:
            n_t (int): number time steps
            n_m (int): number realizations
            dt (float): time step size
            U (float): mean wind speed
            sig (float): wind standard deviation
            tau (float): Kaimal integral time scale
            rho (float): concentration parameter
            mu (float): location parameter
        
        Returns:
            t (numpy array): time
            x (numpy array): [n_t x n_m] array of turbulent wind records
    """
    import numpy as np
    import misc
    import IEC

    scale = 1;    # flag to scale for time discretizaion

    # intermediate parameters
    T    = n_t*dt;                                  # total time
    t    = np.arange(0,n_t)*dt;                     # time vector
    n_f  = misc.uniqueComponents(n_t);              # no. unique components
    df   = 1./T;                                    # frequency resolution
    f    = df*np.arange(0,n_f);                     # frequency array

    if (sig<0):
        print 'Turbulence cannot be less than zero.';
        return []
        
    elif (sig==0):
        x = np.ones(n_t,n_m)*U;
        return [t,x]
        
    else:
        S = IEC.KaimalSpectrum(f,tau,sig);              # Kaimal spectrum
        Sk = S*df;                                  # discrete spectrum
        Xmag = misc.Sk2Xuniq(Sk).reshape(n_f,1)          # S[k] -> |X[k]|
        
        if scale:                                   # spectral scaling
            alpha = misc.spectralScale(Sk,sig,n_t);      #   for discretizaion
            Xmag = alpha*Xmag

        dtheta = wrappedCauchySample(n_f-1,\
            n_m,rho,mu)                             # phase differences
        phases = np.cumsum(dtheta,axis=0)           # phases from f1 to end
        phases = np.append(np.zeros((1,n_m)), \
            phases, axis=0)                         # phases from f0 to end
                
        Xpha = np.exp(1j*phases)                    # phases -> phasors
        
        Xuniq = np.multiply(Xmag,Xpha)              # X = mags*phasors
        Xuniq[0][:] = U;                            # real signal, correct mean
        X = misc.Xuniq2X(Xuniq,n_t)                      # full Fourier vector
        
        x = np.real(np.fft.ifft(X,axis=0) \
            *n_t).reshape(n_t,n_m)                  # time series
                
        return (t,x)