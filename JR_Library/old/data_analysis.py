""" Collection of functions used in data analysis.

Jenni Rinker, 27-April-2015.
"""

def loadmetadata(fname):
    """ Load data and fields from metadata file
    
        Args:
            fname (string): path to file
            
        Returns:
            fields (list): fields in each column
            metadata (numpy array): values for each field and each record
    """
    import numpy as np
    
    # strip fields from header line
    with open(fname,'r') as f:
        fields = f.readline().lstrip('# ').rstrip('\n').split(',')
        
    # read in data table
    metadata = np.loadtxt(fname,delimiter=',',skiprows=1)
    
    return (fields, metadata)
    
def screenmetadata(fields,metadata,dataset):
    """ Screen the metadata for data quality
    
        Args:
            fname (string): name of metadata text file
            dataset (string): toggle for dataset choice
            
        Returns:
            cleandata (numpy array): numpy array with screened data
    """
    import numpy as np
    
    # load unscreened data
#    fields, metadata = loadmetadata(fname)
    
    if (dataset == 'NREL'):
        CSLim  = 3                          # lower cup speed limit
        dir1   = 240                        # CCW edge for direction range
        dir2   = 315                        # CW edge for direction range
        preLim = 2.7                        # lower precipitation limit
        
        # column indices for each value
        CScol  = fields.index('Wind_Speed_Cup')
        dirCol = fields.index('Wind_Direction')
        preCol = fields.index('Precipitation')
#        CScol  = fields.index('Cup_speed')
#        dirCol = fields.index('Wind_direction')
#        preCol = fields.index('Precip_intens')
        
        # filter out the rows with NaN values
        
        
        # screen remaining data
        cleandata = metadata[np.where(metadata[:,CScol] > CSLim)]
        cleandata = cleandata[np.where(cleandata[:,dirCol] >= dir1)]
        cleandata = cleandata[np.where(cleandata[:,dirCol] <= dir2)]
        cleandata = cleandata[np.where(cleandata[:,preCol] >= preLim)]
        
    else:
        print '***ERROR*** That dataset is not coded yet'
        return
        
    return cleandata

def loadNRELmatlab():
    """ Load the metadata table processed in Matlab

        Returns:
            fields (list): names of fields in coumns
            metadata (numpy array): array of atmospheric params
    """
    import scipy.io as scio
    import numpy as np
    
    # path to matlab-processed metadata table
    matpath = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
        '2015-02-28_temporal coherence in data\\code\\dataStruc.mat'

    # load the structure
    struc = scio.loadmat(matpath)

    # extract the information
    dataStruc = struc['dataStruc'][0,0]

    # metadata table
    metadata = dataStruc[0]

    # oddly-formatted field names
    tmp = np.squeeze(dataStruc[1])
        
    # extract readable fieldnames
    fields = [str(ele[0]) for ele in tmp]

    # replace necessary old fieldnames with new ones for screening
    fields[fields.index('Cup_speed')] = 'Wind_Speed_Cup'
    fields[fields.index('Wind_direction')] = 'Wind_Direction'
    fields[fields.index('Precip_intens')] = 'Precipitation'
    fields[fields.index('rho_u')] = 'Concentration_u'

    return (fields, metadata)
    
##def loadsonic(dataset,recordtime):
##    """ Get the time vector and time history for a given field in a dataset
##    
##        Args:
##            dataset (string): toggle to indicate which dataset to query
##            recordtime (float): time at start of record in UTC
##            field (string): desired time history variable
##            
##        Returns:
##            t (numpy array): vector of time increments
##            u (numpy array): longitudinal velocity
##            v (numpy array): transverse velocity
##            w (numpy array): vertical velocity
##            T (numpy array): sonic temperature
##    """
##
##    if (dataset == 'NREL'):
##
##        # 
##
##    else:
##        print '***ERROR*** That dataset is not coded yet.'

def calculateNSAE(x,dist_name,param):
    """ Calculate and return the normalized sum of absolute error

        Args:
            x (numpy array): sample of variables
            dist_name (string): distribution name
            param (tuple): distribution parameters

        Returns:
            NSAE (float): normalized sum of absolute errors between
                empirical and fit CDFs
    """
    import numpy as np
    import scipy.stats
    import matplotlib.pyplot as plt

    x = np.sort(x)                              # resort x to compare CDFs

    N       = x.size                            # no. of elements in sample
    cdf_emp = (np.arange(N)+1.)/(N+1.)            # empirical CDF
    
    dist = getattr(scipy.stats, dist_name)      # define distribution

    cdf_fitted = dist.cdf(x, *param[:-2], \
                loc=param[-2], scale=param[-1]) # theoretical CDF

    NSAE = np.sum(np.abs(cdf_emp-cdf_fitted))/N

    return NSAE
    
def compositeNSAE(x,dist_name,Q,param,param_par):
    """ Return NSAE for composite distribution
    """
    import scipy.stats
    import numpy as np
    
    x       = np.sort(x)                        # resort x to compare CDFs 
    N       = x.size                            # no. of elements in sample
    cdf_emp = (np.arange(N)+1.)/(N+1.)          # empirical CDF
    
    dist     = getattr(scipy.stats, dist_name)    # define distribution
    dist_par = getattr(scipy.stats, 'genpareto')  # define distribution

    # composite
    if Q < 1.0:
        cdf_fitted = np.empty(x.shape)
        cdf_fitted[:N*Q] = dist.cdf(x[:N*Q], *param[:-2], \
            loc=param[-2], scale=param[-1])
        cdf_fitted[N*Q:] = (1- cdf_fitted[N*Q-1])*dist_par.cdf(x[N*Q:]-x[N*Q-1], \
            *param_par[:-2], loc=param_par[-2], scale=param_par[-1]) \
            + cdf_fitted[N*Q-1]        
        
    # single distribution
    else:
        cdf_fitted = dist.cdf(x, *param[:-2], \
                loc=param[-2], scale=param[-1]) # theoretical CDF
                
    NSAE = np.sum(np.abs(cdf_emp-cdf_fitted))/N
    
    return NSAE, cdf_fitted
    

##def calculateNSAE(x, dist_name, param, Q=1.0, param_par=(0.1,0,1)):
##    """ Calculate and return the normalized sum of absolute error for a
##        distribution. For a single distribution, omit Q and the Pareto
##        parameters. Include them for a composite (main + Generalized
##        Pareto) distribution.
##
##        Args:
##            x (numpy array): sample of variables
##            dist_name (string): distribution name
##            param (tuple): main distribution parameters
##            Q (float): cutofff quantile between main and GP distribution
##            param_par (tuple): pareto distribution parameters
##
##        Returns:
##            NSAE (float): normalized sum of absolute errors between
##                empirical and fit CDFs
##    """
##    import scipy.stats
##    import numpy as np
##    
##    x       = np.sort(x)                        # resort x to compare CDFs 
##    N       = x.size                            # no. of elements in sample
##    cdf_emp = (np.arange(N)+1.)/(N+1.)          # empirical CDF
##    
##    dist     = getattr(scipy.stats, dist_name)    # define distribution
##    dist_par = getattr(scipy.stats, 'genpareto')  # define distribution
##
##    # composite
##    if Q < 1.0:
##        cdf_fitted = np.empty(x.shape)          # initialize cdf array
##        cdf_fitted[:N*Q] = dist.cdf(x[:N*Q], *param[:-2], \
##            loc=param[-2], scale=param[-1])     # main part of cdf
##        cdf_fitted[N*Q:] = (1- cdf_fitted[N*Q-1])*dist_par.cdf(x[N*Q:]-x[N*Q-1], \
##            *param_par[:-2], loc=param_par[-2], scale=param_par[-1]) \
##            + cdf_fitted[N*Q-1]                 # GP part of CDF
##        
##    # single distribution
##    else:
##        cdf_fitted = dist.cdf(x, *param[:-2], \
##                loc=param[-2], scale=param[-1]) # theoretical CDF
##                
##    NSAE = np.sum(np.abs(cdf_emp-cdf_fitted))/N # calculate NSAE
##    
##    return NSAE, cdf_fitted
##
##
##def maincompositeNSAE(x,Q,dist_name,param):
##    """ Get normalized sum of absolute error for main distribution
##    """
##    import scipy.stats
##    import numpy as np
##    
##    # convert param to tuple
##    param = tuple(np.ndarray.tolist(param))
##    
##    x = np.sort(x)                              # resort x to compare CDFs
##
##    N       = x.size                            # no. of elements in sample
##    N_calc  = int(Q*N)                          # no. elements in calculation
##    cdf_emp = (np.arange(N)+1.)/(N+1.)          # empirical CDF
##    
##    dist = getattr(scipy.stats, dist_name)      # define distribution
##    
##    if (dist_name == 'genpareto'):        
##        cdf_fitted = dist.cdf(x, *param[:-2], \
##                    loc=0, scale=param[-1]) # theoretical CDF        
##    else:
##        cdf_fitted = dist.cdf(x, *param[:-2], \
##                    loc=param[-2], scale=param[-1]) # theoretical CDF
##                
##    NSAE = np.sum(np.abs(cdf_emp[:N_calc]-cdf_fitted[:N_calc]))/N_calc
##
##    return NSAE
##
##        
##def fitcompositedistribution(x,Q,dist_name):
##    """ Return a tuple of distribution parameters (organized into two tuples,
##        first is for main distribution and second is for gen. pareto).
##        
##        If Q is zero, returns empty tuple for Gen Pareto
##    """
##    import scipy.stats
##    import numpy as np
##    from scipy.optimize import minimize
##    
##    # fit data, initialize parameters
##    dist   = getattr(scipy.stats, dist_name)
##    param0 = dist.fit(x)
##    
##    # function to optimize for main distribution
##    fun = lambda param: maincompositeNSAE(x,Q,dist_name,param)
##    
##    # optimize main distribution parameters
##    res = minimize(fun,param0)
##    param_opt = tuple(np.ndarray.tolist(res.x))
##    
##    if not res.success:
##        print Q, dist_name, res.message
##    
##    # if we're fitting a composite distribution
##    if Q < 1.0:
##        
##        # fit generalized Pareto if Q < 1.0
##        x = np.sort(x)
##        N = x.size
##        x_par = x[np.where(x > x[N*Q])] - x[N*Q]
##        
##        # fit data, initialize parameters
##        dist_par   = getattr(scipy.stats, 'genpareto')
##        param0_par = dist_par.fit(x_par,loc=0)
##        
##        # function to optimize for GP distribution
##        fun_par = lambda param: maincompositeNSAE(x_par,1.0,'genpareto',param)
##        
##        # optimize GP parameters
##        res_par = minimize(fun_par,param0_par)
##        param_opt_par = tuple(np.ndarray.tolist(res_par.x))
##        
##        if not res.success:
##            print(Q, dist_name, res.message)
##        
##    # if it's a single distribution
##    else:
##        
##        param_opt_par = ()
##    
##    return (param_opt,param_opt_par)
    
