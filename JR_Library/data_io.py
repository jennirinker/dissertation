""" Subfunctions for loading data files into Python for analyses.

Jenni Rinker, 06-Apr-2015.
"""

def dateToM4fname(datetime):
    """ Convert array of date information to wildacrd name for M4
        dataset.

        Args:
            dateVec (list): [year,month,day,hour,minute], strings

        Returns:
            fname (string): beginning of filename
    """
    
    fname = datetime[1] + '_' + datetime[2] + '_' + datetime[0] \
            + '_' + datetime[3] + '_' + datetime[4]
    
    return fname
