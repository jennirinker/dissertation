""" Subfunctions for loading TurbSim/data files into Python for analyses.

Jenni Rinker, 26-Mar-2015.
"""

class tsin:
    """ A new class to track TurbSim input parameters of interest.
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

def readInput(filename):
    """ Read specific parameters of interest from TurbSim .inp file and save in
        custom Python class tsin
        
        Args:
            filename (string): path to TS file to read
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