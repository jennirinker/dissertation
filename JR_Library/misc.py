""" Collection of miscellaneous functions used in my analyses.

Jenni Rinker, 26-Mar-2015.
"""

def wrap(theta):
    """ Wrap angle to [-pi,pi).
    
        Args:
            theta (float): angle in radians
    
        Returns:
            theta_out (float): wrapped angle
    """
    import numpy as np
    
    theta_out = np.multiply(theta<=np.pi,theta) + \
        np.multiply(theta>np.pi,theta-2*np.pi)  # subtract 2pi from [pi,2pi)
    
    return theta_out

def uniqueComponents(n_t):
    """ Number of unique Fourier components for a vector of length n_t.

        Args:
           n_t (int): vector length
        
        Returns:
           n_f (int): no. unique components
    """
    import numpy as np

    n_f = int(np.ceil((float(n_t)-1) \
        /2)+1);

    return n_f

def X2Sk(X):
    """ Convert Fourier array to discrete, one-sided spectral mass values.
        
        Args:
            X (numpy array): Fourier vector
        
        Returns:
            Sk (numpy array): array of discrete spectral values
    """
    import numpy as np

    n_t = X.size;                           # no. elements
    n_f = uniqueComponents(n_t);            # no. unique comps
    
    Xuniq = X[:n_f];                        # unique Fourier comps
    
    Sk = np.real(2*np.multiply(Xuniq, \
            np.conj(Xuniq)));               # discrete spectrum

    return Sk

def Sk2Xuniq(Sk):
    """ Convert array of discrete spectral values in Sk to array of unique 
        Fourier magnitudes X.

        Args:
            Sk (numpy array): array of discrete spectral values
        
        Returns:
            Xuniq (numpy array): array of unique Fourier components
    """
    import numpy as np

    Xuniq = np.sqrt(Sk/2);

    return Xuniq

def Xuniq2X(Xuniq,n_t):
    """ Convert array of unique Fourier components to full Fourier array.
        
        Args:
            Xuniq (numpy array): array of unique Fourier components
            n_t (int): no. of components of full array
        
        Returns:
            X (numpy array): array of all Fourier components
    """
    import numpy as np

    if (n_t % 2):                           # even no. of elements
        Xuniq[0] = np.abs(Xuniq[0])         # ensure real signal
        X = np.concatenate((Xuniq,np.conj( \
            Xuniq[1:][::-1])),axis=0);
    else:                                   # odd no. of elements
        Xuniq[0] = np.abs(Xuniq[0])         # ensure real signal
        Xuniq[-1] = np.abs(Xuniq[-1])       # ensure real signal
        X = np.concatenate((Xuniq,np.conj( \
                    Xuniq[1:-1][::-1])),axis=0);

    return X

def spectralScale(Sk,sig,n_t):
    """ Spectral scaling factor for discrete power spectral coefficients Sk.
        
        Args:
            Sk (numpy array): array of discrete spectral values
            sig (float): desired standard deviation of time history
            n_t (int): number of components in time history
        
        Returns:
            alpha (float): spectral scaling factor
    """
    import numpy as np
    
    # if (2+)D array is fed in, halt with error
    if ((len(Sk.shape)>1) and (Sk.shape[0] != 1 and Sk.shape[1] != 1)):
        print 'ERROR: spectralScale only works on 1D arrays'
        return []

    Xuniq   = Sk2Xuniq(Sk)                      # unique Fourier components
    Xuniq[0] = 0                                # neglect DC component
    X       = Xuniq2X(Xuniq,n_t)                # entire Fourier vector
    sumXkSq = np.sum(np.power(np.abs(X),2));    # sum |X[k]|^2
    alpha   = np.sqrt((n_t-1.)/n_t \
        *(sig**2)/sumXkSq);                     # scaling factor
    
    return alpha

def wrappedCauchyPDF(theta,rho,mu):
    """ Probability density function of wrapped Cauchy distribution.
    
        Args:
            theta (numpy array): angles (rads) to evaluate PDF
            rho (float): concentration parameter
            mu (float): location parameter
       
        Returns:
            f (numpy array): evaluation of wrapped Cauchy PDF at angles in theta
    """
    import numpy as np
    
    if ((rho < 0) or (rho>1)):
        print 'rho must be between 0 and 1'
        return []
    
    f = 1/(2*np.pi)*(1 - rho**2) / \
        (1 + (rho**2) - 2*rho*np.cos(theta-mu))
    
    return f