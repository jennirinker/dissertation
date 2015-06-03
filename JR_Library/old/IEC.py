""" Collection of functions to return values from IEC 61400-1 Edition 3.

Jenni Rinker, 26-Mar-2015.
"""

def VelProfile(z,zhub,Vhub):
    """ IEC velocity profile for values in array z with hub height zhub and hub
        velocity Vhub.
        
        Args:
            z (numpy array): array of heights in meters for velocity calcs
            zhub (float): hub-height of turbine in meters
            Vhub (float): hub-height mean velocity in m/s
        
        Returns:
            V (numpy array): array of mean wind speeds at heights in z
    """
    import numpy as np

    alpha = 0.2;

    V = Vhub * np.power( \
        z/zhub, alpha );

    return V

def IrefFromClass(turbc):
    """ Turbulence intensity for turbulence class turbc
    
        Args:
            turbc (string): turbulence class from IEC 61400-1 Ed. 3
        
        Returns:
            Iref (float): reference turbulence intensity
    """

    if ('A' in turbc):                  # class A
        Iref = 0.16;
    elif ('B' in turbc):                # class B
        Iref = 0.14;
    elif ('C' in turbc):                # class C
        Iref = 0.12;
    else:                               # user's choice
        Iref = float(turbc);

    return Iref

def calcSigma1(Iref,Vhub):
    """ Longitudinal standard deviation
    
        Args:
            Iref (float): reference hub-height turbulence intensity
            
        Returns:
            sigma1 (float): longitudinal standard deviation
    """
    
    sigma1 = Iref*(0.75 * Vhub + 5.6);
    
    return sigma1

def TiProfile(z,zhub,Vhub,turbc):
    """ Turbulence intensity profile at heights z for hub height zhub, 
        hub-height wind speed Vhub, and turbulence class turbc.
        
        Args:
            z (numpy array): array of heights in meters for velocity calcs
            zhub (float): hub-height of turbine in meters
            Vhub (float): hub-height mean velocity in m/s
            turbc (string): turbulence class from IEC 61400-1 Ed. 3
            
        Returns:
            Ti (numpy array): turbulence intensity at heights in z
    """

    
    Iref = IrefFromClass(turbc);                # reference hub-height TI  
    sigma1 = calcSigma1(Iref,Vhub)              # reference standard deviation
    V = VelProfile(z,zhub,Vhub);                # mean wind profile
    Ti = sigma1 / V;                            # TI profile
    
    return Ti

def calcLambda1(zhub):
    """ Longitudinal scale parameter
    
        Args:
            zhub (float): hub-height of turbine in meters
        
        Returns:
            Lambda1 (float): longitudinal scale parameter
    """
    
    Lambda1 = (0.7 * zhub)*(zhub <= 60) \
        + 42*(zhub > 60);               # longitudinal scale parameter
    
    return Lambda1

def SpatialCoherence(zhub,Vhub,rsep,f):
    """ Spatial coherence function for hub-height wind speed Vhub and hub 
        height zhub at separation distance rsep and frequencies f.
        
        Args:
            zhub (float): hub-height of turbine in meters
            Vhub (float): hub-height mean velocity in m/s
            rsep (float): separation distance in m
            f (numpy array): frequencies in Hz
        
        Returns:
            Coh (numpy array): values of spatial coherence function
    """
    import numpy as np

    Lambda1 = calcLambda1(zhub)         # longitudinal scale parameter
    Lc = 8.1*Lambda1;                   # coherence scale parameter

    Coh = np.exp(-12*np.sqrt( \
        np.power(f*rsep/Vhub,2) + \
        np.power(0.12*rsep/Lc,2) ));    # spatial coherence

    return Coh

def KaimalSpectrum(f,tau,sig):
    """ Kaimal spectrum (continuous, one-sided) for frequency f and time
        length scale tau = L/U.

        Args:
            f (numpy array): frequencies
            tau (float/numpy array): integral time scale (L/U)
            sig (float): standard deviation
       
        Returns:
            S (numpy array): Kaimal spectrum evaluated at f, tau, sig
    """
    import numpy as np

    S = (sig**2)*(4.*tau)/ \
        np.power(1.+6.*f*tau,5./3.);            # Kaimal 1972

    return S

def PSDs(zhub,Vhub,turbc,f):
    """ Continuous power spectral densities for all 3 turbulence components
    
        Args:
            zhub (float): hub-height of turbine in meters
            Vhub (float): hub-height mean velocity in m/s
            turbc (string): turbulence class from IEC 61400-1 Ed. 3
            f (numpy array): frequencies in Hz
            
        Returns:
            Su (numpy array): longitudinal PSD
            Sv (numpy array): lateral PSD
            Sw (numpy array): vertical PSD
    """

    Iref = IrefFromClass(turbc);        # referene turbulence intensity
    sigma1 = calcSigma1(Iref,Vhub);     # longitudinal standard deviation
    sigma2 = 0.8*sigma1;                # lateral standard deviaiton
    sigma3 = 0.5*sigma1;                # vertical standard deviation
    Lambda1 = calcLambda1(zhub);        # longitudinal scale parameter
    L1 = 8.1*Lambda1;                   # longitudinal integral scale   
    L2 = 2.7*Lambda1;                   # lateral integral scale  
    L3 = 0.66*Lambda1;                  # vertical integral scale

    Su = KaimalSpectrum(f,L1/Vhub,sigma1);  # longitudinal spectrum
    Sv = KaimalSpectrum(f,L2/Vhub,sigma2);  # lateral spectrum
    Sw = KaimalSpectrum(f,L3/Vhub,sigma3);  # vertical spectrum

    return (Su,Sv,Sw)
