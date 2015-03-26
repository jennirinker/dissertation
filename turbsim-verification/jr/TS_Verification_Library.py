
class tsin:
    """ A new class to track TurbSim input parameters
             of interest.
    """
    
    def __init__(self):
        self.ScaleIEC = [];
        self.n_z      = [];
        self.n_y      = [];
        self.dt       = [];
        self.t_anal   = [];
        self.t_use    = [];
        self.zhub     = [];
        self.turbc    = [];
        self.zref     = [];
        self.Vref     = [];

def readInput( filename ):
    """ Read specific parameters of interest from
             TurbSim .inp file and save in custom Python
             class `tsin
    """

    if ('.inp' in filename):

        # define line numbers for each parameter
        line_ScaleIEC = 15;
        line_n_z      = 18;
        line_n_y      = 19;
        line_dt       = 20;
        line_t_anal   = 21;
        line_t_use    = 22;
        line_zhub     = 23;
        line_turbc    = 32;
        line_zref     = 36;
        line_Vref     = 37;

        # initialize the class variable
        tsinput = tsin();

        # open file and extract parameters
        with open( filename, 'r' ) as f:
            idx = 1;
            for line in f:
                if ( idx == line_ScaleIEC ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.ScaleIEC = int(val);
                elif ( idx == line_n_z ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.n_z = int(val);
                elif ( idx == line_n_y ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.n_y = int(val);
                elif ( idx == line_dt ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.dt = float(val);
                elif ( idx == line_t_anal ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.t_anal = float(val);
                elif ( idx == line_t_use ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.t_use = float(val);
                elif ( idx == line_zhub ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.zhub = float(val);
                elif ( idx == line_turbc ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.turbc = val;
                elif ( idx == line_zref ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.zref = float(val);
                elif ( idx == line_Vref ):
                    val = line.lstrip(' ').split(' ')[0]\
                          .replace('\"','');
                    tsinput.Vref = float(val);
                idx += 1;

        return tsinput;

    elif ('.sum' in filename):

        print 'Cannot interpret .sum yet';
        return [];

    else:

        print 'That is neither a .sum or .inp';
        return [];

def calculateTurbSimSC(fname,rsep):
    """ Loads turbsim file in fname and returns
        the spatial coherence at frequency f
        given separation distance rsep.
    """

    import pyts.io.main as io
    import matplotlib.pyplot as plt
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
    [f,Cohij] = spatialCoherence(Xi_all,Xj_all,df);

    return [f, Cohij]

def spatialCoherence(Xi,Xj,df):
    """ Takes [n_f] x [n_p] arrays Xi and Xj, where
        n_p is the number of pairs of points to be
        used in the calculation and returns arrays
        of the frequency and coherence.
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

    return [f,Cohij]

def velocityProfile(tsout):
    """ An array of the mean velocities evaluated
    at each point in the grid.
    """
    import numpy as np

    # rename grid for simplicity
    grid = tsout.grid;

    # initialize U array
    U = np.empty( grid.shape );

    # loop through grid pts, evaluate U, and save
    for iy in range(0,grid.n_y):
	for iz in range(0,grid.n_z):
		U[iz,iy]  = np.mean( tsout.u[iz,iy,:] );

    return U

def hubHeightPSDs(tsout):
    """ Frequency vector and hub-height spectra
        for TurbSim output.
    """
    import numpy as np

    # rename grid for simplicity
    grid = tsout.grid;  

    n_t = tsout.uhub.size;          # stored grid.n_t incorrect
    n_f = np.ceil((n_t-1)/2+1);     # unique Fourier (counts DC)

    # Fourier transform
    Xu = np.fft.fft( tsout.uhub )[0:n_f]/n_t;
    Xv = np.fft.fft( tsout.vhub )[0:n_f]/n_t;
    Xw = np.fft.fft( tsout.whub )[0:n_f]/n_t;

    # Power spectral densities
    Su = np.real( np.multiply( Xu, 
            np.conj( Xu ) ) );
    Sv = np.real( np.multiply( Xv, 
            np.conj( Xv ) ) );
    Sw = np.real( np.multiply( Xw, 
            np.conj( Xw ) ) );

    return [Su,Sv,Sw]

def hubHeightPhaseDiffs(tsout):
    """ Hub-height phase differences.
    """
    import numpy as np

    # rename grid for simplicity
    grid = tsout.grid;  

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

    return [dthetau,dthetav,dthetaw]

def IEC_VelProfile(z,tsout):
    """ IEC velocity profile for TurbSim 
        file tsout at heights in z.
    """
    import numpy as np

    alpha = 0.2;            # IEC 61400-1 Ed-3
    ihubz = tsout.ihub[0];  # hub height

    V = tsout.UHUB * np.power( \
        z/tsout.z[ihubz], alpha );

    return V

def IEC_Iref(tsin):
    """ Return turbulence intensity for
    TurbSim input info in tsin.
    """

    if ('A' in tsin.turbc):     # class A
        Iref = 0.16;
    elif ('B' in tsin.turbc):   # class B
        Iref = 0.14;
    elif ('C' in tsin.turbc):   # class C
        Iref = 0.12;
    else:                       # user's choice
        Iref = float(tsin.turbc);

    return Iref

def IEC_TiProfile(z,tsin,tsout):
    """ IEC turbulence intensity profile
        for TurbSim file tsout at heights
        in z.
    """

    # expected hub-height TI for 15 m/s
    Iref = IEC_Iref(tsin);

    # reference standard deviation
    sigma1 = Iref*(0.75 * \
        tsout.UHUB + 5.6);

    # mean wind profile
    V = IEC_VelProfile(z,tsout);

    # TI profile
    Ti = sigma1 / V;
    
    return Ti

def IEC_SpatialCoherence(rsep,f,tsout):
    """ IEC spatial coherence function
        for TurbSim file tsout at separation
        distance rsep and frequencies f.
    """
    import numpy as np

    # intermediate parameters
    Vhub  = tsout.UHUB;     # hub-height wind speed
    ihubz = tsout.ihub[0];  # hub-height index
    zhub = tsout.z[ihubz];  # hub height
    Lambda1 = (0.7*zhub)*(zhub<=60) \
        + 42*(zhub>60);     # longitudinal scale parameter
    Lc = 8.1*Lambda1;       # coherence scale parameter

    Coh = np.exp(-12*np.sqrt( \
        np.power(f*rsep/Vhub,2) + \
        np.power(0.12*rsep/Lc,2) ));    # spatial coherence

    return Coh
    
def IEC_PSDs(f,tsin,tsout):
    """ Spectra for IEC model.
    """
    import numpy as np

    # intermediate parameters
    df = tsout.grid.df;
    Vhub  = tsout.UHUB;     # hub-height wind speed
    ihubz = tsout.ihub[0];  # hub-height index
    zhub = tsout.z[ihubz];  # hub height
    Iref = IEC_Iref(tsin);  # reference TI
    sigma1 = Iref * (0.75 * \
        Vhub+ 5.6);
    sigma2 = 0.8*sigma1;
    sigma3 = 0.5*sigma1;
    Lambda1 = (0.7*zhub)*(zhub<=60) \
        + 42*(zhub>60);     # longitudinal scale parameter
    L1 = 8.1*Lambda1;       # longitudinal integral scale   
    L2 = 2.7*Lambda1;       # lateral integral scale  
    L3 = 0.66*Lambda1;      # vertical integral scale  

    Su = (sigma1**2)*(4*L1/Vhub)/ \
         np.power(1+6*f*L1/Vhub,5./3.) * df;
    Sv = (sigma2**2)*(4*L2/Vhub)/ \
         np.power(1+6*f*L2/Vhub,5./3.) * df;
    Sw = (sigma3**2)*(4*L3/Vhub)/ \
         np.power(1+6*f*L3/Vhub,5./3.) * df;

    return [Su,Sv,Sw]
