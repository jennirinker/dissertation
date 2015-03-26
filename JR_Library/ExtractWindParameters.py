""" Collection of functions to extract parameters from time series.

Jenni Rinker, 26-Mar-2015.
"""

def samplePhaseCoherence(theta):
    """ Return concentration and location parameters for a sample of wrapping
        random variables.
        
        Args:
            theta (numpy array): 1D array of sample of angles
       
        Returns:
            rho (float): concentration parameter
            mu (float): location parameter
    """
    import numpy as np
    
    # if (2+)D array is fed in, halt with error
    if ((len(theta.shape)>1) and (theta.shape[0] != 1 and theta.shape[1] != 1)):
        print 'ERROR: samplePhaseCoherence only works on 1D arrays'
        return []
    
    n_t = theta.size                        # number of elements in sample
    V = np.sum(np.exp(1j*theta))/n_t        # mean resultant vector
    rho = np.abs(V)                         # concentration parameter
    mu  = np.angle(V)                       # location parameter
    
    return (rho,mu)

def signalPhaseCoherence(x):
    """ Return concentration and location parameters for a time history.
    
        Args:
            x (numpy array): 1D array of time history
       
        Returns:
            rho (float): concentration parameter
            mu (float): location parameter
    """
    import numpy as np
    import misc
    
    # if (2+)D array is fed in, halt with error
    if ((len(x.shape)>1) and (x.shape[0] != 1 and x.shape[1] != 1)):
        print 'ERROR: signalPhaseCoherence only works on 1D arrays'
        return []
    
    n_t = x.shape[0]                        # no. of total components
    n_f = misc.uniqueComponents(n_t)        # no. of unique components
    
    X = np.fft.fft(x,axis=0)/n_t            # Fourier vector
    Xuniq = X[:n_f]                         # unique Fourier components
    dtheta = np.angle(np.divide( \
        Xuniq[1:],Xuniq[:-1]))              # phase differences
    rho, mu = samplePhaseCoherence( \
        dtheta)                             # temporal coherence parameters
    
    return (rho,mu)

def signal2Sk(x):
    """ Return array of discrete spectrum for signal x.

        Args:
            x (numpy array): array of time hietory
       
        Returns:
            Sk (numpy array): array of discrete spectral values
    """
    import numpy as np
    import misc

    n_t  = x.shape[0];                      # no. of elements
    X = np.fft.fft(x,axis=0)/n_t;           # Fourier vector
    Sk = misc.X2Sk(X);                      # discrete PSD

    return Sk

def spatialCoherence(Xi,Xj,df):
    """ Takes [n_f] x [n_p] arrays Xi and Xj, where n_p is the number of pairs 
        of points to be used in the calculation and returns arrays of the 
        frequency and coherence.
        
        Args:
            Xi ([n_f x n_p] numpy array): Fourier coefficients at first point
            Xj ([n_f x n_p] numpy array): Fourier coefficients at second point
    """

    import numpy as np

    # generate frequency vector
    n_f = Xi.shape[0];
    f   = np.arange(0,n_f)*df;

    # calculate expected values
    EXiXi = np.real(np.mean(np.multiply(Xi,\
        np.conj(Xi)),axis=1));
    EXjXj = np.real(np.mean(np.multiply(Xj,\
        np.conj(Xj)),axis=1));
    EXiXj = np.mean(np.multiply(Xi,\
        np.conj(Xj)),axis=1);

    # calculate PSDs/CSDs
    Sii = 2 * EXiXi / df;
    Sjj = 2 * EXjXj / df;
    Sij = 2 * EXiXj / df;

    # calculate coherence
    Cohij = np.abs( np.divide( Sij, \
        np.multiply(np.sqrt(Sii),np.sqrt(Sjj)) ) );

    return (f, Cohij)

def calculateKaimal(x,dt):
    """ Return the optimal Kaimal time scale tau found using a grid search. 
        Can be converted to length scale using L = tau * U.

        Args:
            x (numpy array): time history
            dt (float): time step
        
        Returns:
            tau (float): optimal Kaimal length scale
    """
    import numpy as np
    import misc
    import IEC
    
    # if (2+)D array is fed in, halt with error
    if ((len(x.shape)>1) and (x.shape[0] != 1 and x.shape[1] != 1)):
        print 'ERROR: calculateKaimal only works on 1D arrays'
        return []

    # grid search parameters
    n_m = 200;                                  # number of points in grid
    tau_l = -1;                                 # left tau coefficient
    tau_r = 3;                                  # right tau coefficient

    # intermediate parameters
    n_t  = x.shape[0]                           # no. of time steps
    T    = n_t*dt                               # total time
    n_f  = misc.uniqueComponents(n_t)           # no. unique components
    df   = 1/T                                  # frequency resolution
    f    = df*np.arange(0,n_f). \
        reshape(n_f,1)                          # frequency vector
    sig  = np.std(x,ddof=1)                     # std. deviation

    X_dat  = np.fft.fft(x,axis=0). \
        reshape(n_t,1)/n_t                      # Fourier vector of data
    Sk_dat = misc.X2Sk(X_dat)[1:]* \
        np.ones((1,n_m))                        # data discrete PSD from f1 up

    taugrid = np.logspace(tau_l,tau_r,n_m). \
        reshape(1,n_m)                          # grid of tau for search
    
    S_kaim = IEC.KaimalSpectrum(f,taugrid,sig)  # continuous Kaimal spectrum
    Sk_kaim = S_kaim * df                       # discrete Kaimal spectrum
    Sk_theo = np.empty(Sk_kaim.shape)
    for i in range(n_m):
        alpha = misc.spectralScale(Sk_kaim[:,i],sig,n_t)
        Sk_theo[:,i] = (alpha**2)*Sk_kaim[:,i]
    Sk_theo = Sk_theo[1:]                       # strip off DC component
    
    J = np.sum(np.power(Sk_dat-Sk_theo, \
        2),axis=0)                              # array of squared errors

    i_min = np.argmin(J)                        # index of min error

    tau = taugrid[0,i_min]                      # tau_minError
    
    return tau

def calculateTurbSimSC(fname,rsep):
    """ Loads turbsim file in fname and returns the spatial coherence at 
        frequency f given separation distance rsep.
        
        Args:
            fname (string): filename
            rsep (float): separation distance in meters
            
        Returns:
            
    """
    import pyts.io.main as io
    import numpy as np

    # read file
    tsout = io.readModel( fname );

    # extract grid
    grid = tsout.grid;

    # define useful parameters for simpler code
    n_p = grid.n_p;
    n_t = tsout.uhub.size;          # grid.n_t is not correct
    n_f = np.ceil((n_t-1)/2+1);     # unique components counting DC
    df  = grid.df;

    # initialize PSD/CSD arrays
    Xi_all = np.zeros( [n_f, 1]);
    Xj_all = np.zeros( [n_f, 1]);

    # loop through points in grid
    for ii in range(0,n_p):
        for jj in range(ii,n_p):

            # check if correct separation
            if ( grid.dist(ii,jj) == rsep ):

                # get grid indices
                [iiz,iiy] = grid.ind2sub(ii);
                [jjz,jjy] = grid.ind2sub(jj);

                # extract time histories
                ui = tsout.u[iiz,iiy,:];
                uj = tsout.u[jjz,jjy,:];

                # calculate Fourier coefficients
                Xi = np.fft.fft(ui)[0:n_f].reshape([n_f,1])/n_t;
                Xj = np.fft.fft(uj)[0:n_f].reshape([n_f,1])/n_t;

                # append to array
                Xi_all = np.hstack([Xi_all, Xi]);
                Xj_all = np.hstack([Xj_all, Xj]);

    # remove zeros from initial column
    Xi_all = Xi_all[:,1:];
    Xj_all = Xj_all[:,1:];

    # calculate spatial coherence of arrays
    f, Cohij = spatialCoherence(Xi_all,Xj_all,df);

    return (f, Cohij)

def TurbSimVelProfile(fname):
    """ An array of the mean velocities evaluated at each point in the TurbSim
        grid of file fname.
        
        Args:
            fname (string): TurbSim filename
            
        Returns:
            U (numpy array): mean velocities at heights in grid
    """
    import pyts.io.main as io
    import numpy as np

    # read file
    tsout = io.readModel(fname);

    # rename grid for simplicity
    grid = tsout.grid;

    # initialize U array
    U = np.empty( grid.shape );

    # loop through grid pts, evaluate U, and save
    for iy in range(0,grid.n_y):
	for iz in range(0,grid.n_z):
		U[iz,iy]  = np.mean( tsout.u[iz,iy,:] );

    return U

def TurbSimHHPSDs(fname):
    """ Frequency vector and hub-height spectra for TurbSim output.
    
        Args:
            fname (string): TurbSim filename
            
        Returns:
            Suk (numpy array): discrete longitudinal PSD
            Svk (numpy array): discrete lateral PSD
            Swk (numpy array): discrete vertical PSD
    """
    import pyts.io.main as io

    # read file
    tsout = io.readModel(fname);

    # calculate PSDs
    Suk = signal2Sk(tsout.uhub)
    Svk = signal2Sk(tsout.vhub)
    Swk = signal2Sk(tsout.whub)
    
    return (Suk,Svk,Swk)
    
def TurbSimHHPDDs(fname):
    """ Hub-height phase differences
    """
    import pyts.io.main as io
    import numpy as np

    # read file
    tsout = io.readModel(fname);

    n_t = tsout.uhub.size;          # stored grid.n_t incorrect
    n_f = np.ceil((n_t-1)/2+1);     # unique Fourier (counts DC)

    # Fourier transform
    Xu = np.fft.fft( tsout.uhub )[0:n_f]/n_t;
    Xv = np.fft.fft( tsout.vhub )[0:n_f]/n_t;
    Xw = np.fft.fft( tsout.whub )[0:n_f]/n_t;

    # Phase differences
    dthetau = np.angle( np.divide(Xu[1:], Xu[:-1]) );
    dthetav = np.angle( np.divide(Xv[1:], Xv[:-1]) );
    dthetaw = np.angle( np.divide(Xw[1:], Xw[:-1]) );

    return (dthetau,dthetav,dthetaw)