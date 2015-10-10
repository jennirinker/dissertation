"""
Temp testing reading text files copied from excel
"""


def mygenfromtxt(fname,header=True,units=False,delimiter='\t'):
    """ Load delimited array into numpy with optional header and units. Assumes
        units line will only be present if header line present.
    
        Args:
            fname (string): filename or path to file
            header (Boolean): whether to read header line
            units (Boolean): whether to read units line
            delimiter (string): string used to separate values
            
        Returns:
    """
    import numpy as np
    
    # read header and unit lines if appropriate
    with open(fname,'r') as f:
        
        # read header line if present
        if header:
            header_text = f.readline().strip('\n').split('\t')
            
            # read units line if present
            if units:
                units_text = f.readline().strip('\n').split('\t')
            else:
                units_text = []
                
        else:
            header_text = []
            
    # use numpy to load numeric data
    data = np.genfromtxt(fname,delimiter=delimiter,
                         skip_header=header+units)
    
    # return different values based on inputs
    if not header:
        return data
    elif (header and (not units)):
        return data, header_text
    else:
        return data, header_text, units_text


fname = 'textfiles//0.75_Blade.txt'
header = 1
units  = 0

data, header = mygenfromtxt(fname)