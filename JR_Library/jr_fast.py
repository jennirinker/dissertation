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
    
def CreateFAST7Dict(fast_fpath):
    """ Build and save FAST 7 Python dictionary from input file
    
        Args:
            fast_fpath (string): path to .fst file
            
        Returns:
            turb_dict (dictionary): dictionary of turbine parameters
    """
    
    # ensure path is to a .fst file 
    if not fast_fpath.endswith('.fst'):
        err_str = 'Path {:s} is not to a FAST 7 input file'.format(fast_fpath)
        ValueError(err_str)
        
    # change directories to turbine directory
    turb_dir = os.path.dirname(fast_fpath)
    os.chdir(turb_dir)
    
    # ====================== initialize dictionary ============================
    turb_dict = {}
    turb_dict['TurbName'] = os.path.basename(fast_fpath)
    turb_dict['TurbDir']  = turb_dir
    
    # ==================== read data from .fst file ===========================
    with open(fast_fpath,'r') as f:
        
        # read first four lines manually
        f.readline()
        f.readline()
        turb_dict['FstCmt1'] = f.readline().rstrip('\n')
        turb_dict['FstCmt2'] = f.readline().rstrip('\n')
        
        # read up to OutList automatically
        value = ''
        while ( value != 'OutList'):
            line = f.readline()
            
            # if line doesn't start with dashes, it is a parameter
            if ( line[:2] != '--' ):
                # convert to float if number
                try:
                    value = float(line.split()[0])
                # otherwise remove quotes if present
                except ValueError:
                    value = line.split()[0].rstrip('\"').lstrip('\"')
                key   = line.split()[1]
                turb_dict[key] = value

        # read OutList automatically but save differently than above
        OutList = []
        line = f.readline()
        while ( line[:3] != 'END' ):
            OutList.append(line)
            line = f.readline()
        turb_dict['OutList'] = OutList
        
    # ============== read data from platform file if used =====================
    if turb_dict['PtfmModel']:
        with open(turb_dict['PtfmFile'],'r') as f:
            
            # read first four lines manually
            f.readline()
            f.readline()
            line = f.readline().rstrip('\n')
            turb_dict['PtfmCmt'] = line
            
            # read remaining lines automatically
            line = f.readline().rstrip('\n')
            while ( line ):
                
                # if line doesn't start with dashes, it's a parameter
                if ( line[:2] != '--' ):
                    # convert to float if number
                    try:
                        value = float(line.split()[0])
                    # remove quotes if present
                    except ValueError:
                        value = line.split()[0].rstrip('\"').lstrip('\"')
                    key   = line.split()[1]
                    turb_dict[key] = value
                
                line = f.readline().rstrip('\n')
         
    # =================== read data from tower file ===========================
    with open(turb_dict['TwrFile'],'r') as f:
        
        # read first four lines manually
        f.readline()
        f.readline()
        line = f.readline().rstrip('\n')
        turb_dict['TwrCmt'] = line
        
        # read to HtFract automatically
        value = ''
        while ( value != 'HtFract'):
            line = f.readline()
            
            # if line doesn't start with dashes, it is a parameter
            if ( line[:2] != '--' ):
                # convert to float if number
                try:
                    value = float(line.split()[0])
                # otherwise remove quotes if present
                except ValueError:
                    value = line.split()[0].rstrip('\"').lstrip('\"')
                key   = line.split()[1]
                turb_dict[key] = value
        f.readline()
        
        # read distributed tower properties
        twr_prop = []
        for i_st in range(int(turb_dict['NTwInpSt'])):
            line = f.readline()
            twr_prop.append([float(s) for s in line.rstrip('\n').split()])
        turb_dict['TwrDistProp'] = twr_prop
        
        # read remaining lines automatically
        line = f.readline().rstrip('\n')
        while ( line ):
            
            # if line doesn't start with dashes, it's a parameter
            if ( line[:2] != '--' ):
                # convert to float if number
                try:
                    value = float(line.split()[0])
                # remove quotes if present
                except ValueError:
                    value = line.split()[0].rstrip('\"').lstrip('\"')
                key   = line.split()[1]
                turb_dict[key] = value
            
            line = f.readline().rstrip('\n')  
            
        
    # =============== read data from furling file if used =====================
    if ( turb_dict['Furling'] == 'True' ):
        with open(turb_dict['FurlFile'],'r') as f:
            
            # read first four lines manually
            f.readline()
            f.readline()
            line = f.readline().rstrip('\n')
            turb_dict['FurlCmt'] = line
            
            # read remaining lines automatically
            line = f.readline().rstrip('\n')
            while ( line ):
                
                # if line doesn't start with dashes, it's a parameter
                if ( line[:2] != '--' ):
                    # convert to float if number
                    try:
                        value = float(line.split()[0])
                    # remove quotes if present
                    except ValueError:
                        value = line.split()[0].rstrip('\"').lstrip('\"')
                    key   = line.split()[1]
                    turb_dict[key] = value
                
                line = f.readline().rstrip('\n')
                
    # =================== read data from blade files ===========================
    for i_bl in range(1,int(turb_dict['NumBl'])+1):
        
        # append blade number to all blade keys to differentiate
        bl_str = str(i_bl)
        
        bl_key = 'BldFile({:d})'.format(i_bl)
        with open(turb_dict[bl_key],'r') as f:
            
            # read first four lines manually
            f.readline()
            f.readline()
            line = f.readline().rstrip('\n')
            turb_dict['BlCmt' + bl_str] = line
            
            # read to BlFract automatically
            value = ''
            while ( value != 'BlFract'):
                line = f.readline()
                
                # if line doesn't start with dashes, it is a parameter
                if ( line[:2] != '--' ):
                    # convert to float if number
                    try:
                        value = float(line.split()[0])
                    # otherwise remove quotes if present
                    except ValueError:
                        value = line.split()[0].rstrip('\"').lstrip('\"')
                    key   = line.split()[1] + bl_str
                    turb_dict[key] = value
            f.readline()
            
            # read distributed tower properties
            bld_prop = []
            for i_st in range(int(turb_dict['NBlInpSt' + bl_str])):
                line = f.readline()
                bld_prop.append([float(s) for s in line.rstrip('\n').split()])
            turb_dict['BldDistProp' + bl_str] = twr_prop
            
            # read remaining lines automatically
            line = f.readline().rstrip('\n')
            while ( line ):
                
                # if line doesn't start with dashes, it's a parameter
                if ( line[:2] != '--' ):
                    # convert to float if number
                    try:
                        value = float(line.split()[0])
                    # remove quotes if present
                    except ValueError:
                        value = line.split()[0].rstrip('\"').lstrip('\"')
                    key   = line.split()[1] + bl_str
                    turb_dict[key] = value
                
                line = f.readline().rstrip('\n') 
                
        
    # =================== read data from AeroDyn file ===========================
    with open(turb_dict['ADFile'],'r') as f:
        
        # read first line manually
        line = f.readline().rstrip('\n')
        turb_dict['ADCmt'] = line
        
        # read to NumFoil automatically
        key = ''
        while ( key != 'NumFoil'):
            line = f.readline()
            
            # if line doesn't start with dashes, it is a parameter
            if ( line[:2] != '--' ):
                # convert to float if number
                try:
                    value = float(line.split()[0])
                # otherwise remove quotes if present
                except ValueError:
                    value = line.split()[0].rstrip('\"').lstrip('\"')
                key   = line.split()[1]
                turb_dict[key] = value
        
        # read foil files
        foil_prop = []
        for i_st in range(int(turb_dict['NumFoil'])):
            line = f.readline().split()[0].lstrip('\"').rstrip('\"\n')
            foil_prop.append(line)
        turb_dict['ADFoils'] = foil_prop
        
        # read number of blade nodes
        line = f.readline()     
        # convert to float if number
        try:
            value = float(line.split()[0])
        # otherwise remove quotes if present
        except ValueError:
            value = line.split()[0].rstrip('\"').lstrip('\"') 
        key   = line.split()[1]
        turb_dict[key] = value
        f.readline()
        
        # read blade nodes
        AD_prop = []
        for i_bl in range(int(turb_dict['BldNodes'])):
            line = f.readline().rstrip('\n').split()
            row = []
            for i_col in range(len(line)):
                try:
                    row.append(float(line[i_col]))
                except ValueError:
                    row.append(line[i_col])
            AD_prop.append(row)
        turb_dict['ADDistProp'] = AD_prop
    
    return turb_dict
    
    # Noise file, linearization file, ADAMS file
