#! /usr/bin/python
"""
Create user spectrum file
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
#libpath = '/home/jrinker/git/dissertation/'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np


# define wind/simulation parameters
fname_spc  = 'Testing.spc'
fname_inp   = 'Testing.inp'
TS = {}
TS['URef']  = 10.
TS['sig_u'] = 2.
TS['L_u']   = 340.2
TS['rho_u'] = 0.2
TS['rho_v'] = 0.2
TS['rho_w'] = 0.2
TS['R1']   = 123456
TS['R2']   = 789012

TS['sig_v'] = 0.8*TS['sig_u']
TS['sig_w'] = 0.5*TS['sig_u']
TS['L_v']   = TS['L_u']/8.1*2.7
TS['L_w']   = TS['L_u']/8.1*0.66

TS['n_z']  = 5
TS['n_y']  = 5
TS['DZ']   = 140
TS['DY']   = 140

TS['mu_u']  = np.pi
TS['mu_v']  = np.pi
TS['mu_w']  = np.pi
TS['T']     = 600.
TS['dt']    = 0.05

TS['ZHub'] = 90
TS['ZRef'] = 90

TS['fpath_spc'] = fname_spc

TS['TurbModel'] = 'USRINP'
TS['ProfType']  = 'PL'

wr_dir = ''

jr.WriteTurbSimInputs(fname_inp,TS,wr_dir,
                       fname_spc=fname_spc)
            



