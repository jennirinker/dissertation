"""
A series of Python functions for the creation and analysis of FAST-related
files.

AUTHOR:  Jenni Rinker, Duke University
CONTACT: jennifer.rinker@duke.edu

===============================================================================

NOTES:
    Wind turbine template directory (FAST 7)
    ----------------------------------------
    Code is structured under the assumption that any wind turbine template 
    directories for FAST 7 have non-wind-dependent files (i.e., blade, tower,
    and control files) in the top level of the directory, wind-dependent 
    templates (i.e., AeroDyn and FAST templates) are in subfolder "templates", 
    and steady-state look-up table is in subfolder "steady_state" with the name 
    "<turb_name>_SS.mat".

===============================================================================
"""
import jr_wind
import os
import scipy.io as scio
import numpy as np


def WriteFAST7InputsAll(tmpl_dir,turb_name,wind_dir,
                      **kwargs):
    """ Write FAST 7 input files for all wind files in directory. Assumes wind 
        turbine template directory has organization specified in module 
        docstring. If output directory is specified, all FAST input files are 
        either written or copied to specified directory.
    
        Args:
            tmpl_dir (string): path to wind turbine template directory
            wind_fpaths (list): list of file paths to wind files
            wr_dir (string): directory to write input files to [opt]
            kwargs (dictionary): keyword arguments to WriteFAST7InputsOne [opt]
            
    """
    
    # possible wind file endings
    wind_ends = ('.bts','.wnd')
            
    # get list of wind files from directory
    wind_fnames = [f for f in os.listdir(wind_dir) if f.endswith(wind_ends)]

    # loop through wind files
    for wind_fname in wind_fnames:
        WriteFAST7InputsOne(tmpl_dir,turb_name,wind_fname,
                   wind_dir=wind_dir,
                   **kwargs)
    
    return
    
def WriteFAST7InputsOne(tmpl_dir,turb_name,wind_fname,
                   BlPitch0=None,RotSpeed0=None,
                   wind_dir=None,fileID=None,t_max=630.,
                   wr_dir=None):
    """ Write FAST 7 input files for one wind file. Assumes wind 
        turbine template directory has organization specified in module 
        docstring. If output directory is specified, all FAST input files are 
        either written or copied to specified directory.
    
        Args:
            tmpl_dir (string): path to turbine FAST directory
            turb_name (string): turbine name
            wind_fname (string): name of wind file
            wind_dir (string): path to directory with wind files
            BlPitch0 (list/numpy array): initial blade pitch angles [opt]
            RotSpee0 (float): initial rotor speed [opt]
            fileID (string): file identifier [opt]
            t_max (float): maximum simulation time [opt]
            
    """

    
    # get initial wind speed
    wind_fpath = os.path.join(wind_dir,wind_fname)
    u0 = jr_wind.GetFirstWind(wind_fpath)
    
    print('Writing FAST files for \"{:s}\" '.format(turb_name) + \
            'with wind file {:s}'.format(wind_fpath))
    
    # set optional values as necessary
    GenDOF = 'True'
    if wind_dir is None:
        wind_dir = os.path.join(tmpl_dir,'Wind')
    if wr_dir is None:
        wr_dir = tmpl_dir
    if fileID is None:
        fAD_name  = wind_fname[:-4] + '_AD.ipt'
        fFST_name = wind_fname[:-4] + '.fst'
    else:
        fAD_name  = turb_name+'_'+fileID+'_AD.ipt'
        fFST_name = turb_name+'_'+fileID+'.fst'
    if BlPitch0 is None:
        mdict = scio.loadmat(os.path.join(tmpl_dir,'steady_state',
                                        turb_name+'_SS.mat'),squeeze_me=True)
        LUT        = mdict['SS']
        saveFields = [str(s).strip() for s in mdict['Fields']]
        BlPitch0 = np.interp(u0,LUT[:,saveFields.index('WindVxi')],
                            LUT[:,saveFields.index('BldPitch1')])*np.ones(3)
    if RotSpeed0 is None:
        mdict = scio.loadmat(os.path.join(tmpl_dir,'steady_state',
                                        turb_name+'_SS.mat'),squeeze_me=True)
        LUT        = mdict['SS']
        saveFields = [str(s).strip() for s in mdict['Fields']]
        RotSpeed0 = np.interp(u0,LUT[:,saveFields.index('WindVxi')],
                            LUT[:,saveFields.index('RotSpeed')])
    if t_max is None:
        t_max = jr_wind.GetLastTime(wind_fpath)

    # create filenames
    fAD_temp  = os.path.join(tmpl_dir,'templates',turb_name+'_AD.ipt')
    fAD_out   = os.path.join(wr_dir,fAD_name)
    fFST_temp = os.path.join(tmpl_dir,'templates',turb_name+'.fst')
    fFST_out  = os.path.join(wr_dir,fFST_name)
    
    # write AeroDyn file
    with open(fAD_temp,'r') as f_temp:
        with open(fAD_out,'w') as f_write:
            i_line = 0
            for line in f_temp:
                if i_line == 9:
                    f_write.write(line.format(wind_fpath))
                else:
                    f_write.write(line)
                i_line += 1
                
    # write FAST file
    with open(fFST_temp,'r') as f_temp:
        with open(fFST_out,'w') as f_write:
            i_line = 0
            for line in f_temp:
                if i_line == 9:
                    f_write.write(line.format(t_max))
                elif i_line == 45:
                    f_write.write(line.format(BlPitch0[0]))
                elif i_line == 46:
                    f_write.write(line.format(BlPitch0[1]))
                elif i_line == 47:
                    f_write.write(line.format(BlPitch0[2]))
                elif i_line == 59:
                    f_write.write(line.format(GenDOF))
                elif i_line == 72:
                    f_write.write(line.format(RotSpeed0))
                elif i_line == 160:
                    f_write.write(line.format(fAD_name))
                else:
                    f_write.write(line)
                i_line += 1
                
    return