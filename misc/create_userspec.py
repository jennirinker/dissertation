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

spec = 1
TS   = 1

# define wind/simulation parameters
fspcname  = 'Testing.spc'
fTSname   = 'Testing.inp'
TS = {}
#TS['Uhub']  = 10.
#TS['sig_u'] = 2.
#TS['L_u']   = 340.2
#TS['rho_u'] = 0.2
#TS['rho_v'] = 0.1
#TS['rho_w'] = 0.0
#TS['R1']   = 123456
#TS['R2']   = 789012

#TS['sig_v'] = 0.8*TS['sig_u']
#TS['sig_w'] = 0.5*TS['sig_u']
#TS['L_v']   = TS['L_u']/8.1*2.7
#TS['L_w']   = TS['L_u']/8.1*0.66

TS['n_z']  = 5
TS['n_y']  = 5
TS['DZ']   = 140
TS['DY']   = 140

#TS['mu_u']  = np.pi
#TS['mu_v']  = np.pi
#TS['mu_w']  = np.pi
TS['T']     = 600.
TS['dt']    = 0.05

TS['ZHub'] = 90


    # template filename
    spctemp = 'Template_UsrSpc.spc'
    TStemp  = 'Template_TurbSim.inp'
    
    # convert items in dictionary to variables for coding convenience
    for key in TS.keys():
        locals()[key] = TS[key]
    
    # create unscaled PSDs
    fo, ff, df = 0, 20, 0.0005                          # high-frequency spectra
    fs  = np.arange(fo,ff,df)                           # high-frequency vector
    Su    = jr.KaimalSpectrum(fs,L_u/Uhub,sig_u)        # unscaled u-PSD
    Sv    = jr.KaimalSpectrum(fs,L_v/Uhub,sig_v)        # unscaled v-PSD
    Sw    = jr.KaimalSpectrum(fs,L_w/Uhub,sig_w)        # unscaled w-PSD
    NumF  = fs.size                                     # number of frequencies
    
    # frequenices/spectral values for scaling
    df_s, n_t = 1./T, T/dt                              # simulation parameters
    fs_s = np.arange(jr.uniqueComponents(n_t))*df_s     # simulation freqs
    Su_s  = jr.KaimalSpectrum(fs_s,L_u/Uhub,sig_u)      # simulation u-PSD
    Sv_s  = jr.KaimalSpectrum(fs_s,L_v/Uhub,sig_v)      # simulation v-PSD
    Sw_s  = jr.KaimalSpectrum(fs_s,L_w/Uhub,sig_w)      # simulation w-PSD
    Suk_s, Svk_s, Swk_s = Su_s*df_s,\
                    Sv_s*df_s, Sw_s*df_s                # continuous -> discrete
    alpha1 = jr.spectralScale(Suk_s,sig_u,n_t)**2       # u scale factor
    alpha2 = jr.spectralScale(Svk_s,sig_v,n_t)**2       # v scale factor
    alpha3 = jr.spectralScale(Swk_s,sig_w,n_t)**2       # w scale factor
    Scale1,Scale2,Scale3 = alpha1,alpha2,alpha3         # set scale factors
    
    # create spectral input file
    with open(fspcname,'w') as f_out:
        with open(spctemp,'r') as f_temp:
            
            # print header information
            f_out.write(f_temp.readline())
            f_out.write(f_temp.readline().format(Uhub,sig_u,sig_v,sig_w,
                       L_u,L_v,L_w,rho_u,rho_v,rho_w,mu_u,mu_v,mu_w ))
            f_out.write(f_temp.readline())
            f_out.write(f_temp.readline().format(NumF))
            f_out.write(f_temp.readline().format(Scale1))
            f_out.write(f_temp.readline().format(Scale2))
            f_out.write(f_temp.readline().format(Scale3))
            f_out.write(f_temp.readline())
            f_out.write(f_temp.readline())
            f_out.write(f_temp.readline())
            f_out.write(f_temp.readline())
            
            # print spectra
            for i_f in range(NumF):
                f_out.write('{:>6.4f}        {:>10.6f}        '.format(fs[i_f],Su[i_f]) + \
                        '{:>10.6f}        {:>8.6f}\n'.format(Sv[i_f],Sw[i_f]))
    
    
    # create TurbSim input file
    with open(fTSname,'w') as f_out:
        with open(TStemp,'r') as f_temp: 
    
            i_line = 0
            for line in f_temp:
                if i_line == 4:
                    f_out.write(line.format(R1))
                elif i_line == 5:
                    f_out.write(line.format(R2))
                elif i_line == 18:
                    f_out.write(line.format(n_z))
                elif i_line == 19:
                    f_out.write(line.format(n_y))
                elif i_line == 20:
                    f_out.write(line.format(dt))
                elif i_line == 21:
                    f_out.write(line.format(T))
                elif i_line == 23:
                    f_out.write(line.format(ZHub))
                elif i_line == 24:
                    f_out.write(line.format(DZ))
                elif i_line == 25:
                    f_out.write(line.format(DY))
                elif i_line == 31:
                    f_out.write(line.format(fspcname))
                elif i_line == 38:
                    f_out.write(line.format(ZHub))
                elif i_line == 39:
                    f_out.write(line.format(Uhub))
                elif i_line == 63:
                    f_out.write(line.format(rho_u,mu_u))
                elif i_line == 64:
                    f_out.write(line.format(rho_v,mu_v))
                elif i_line == 65:
                    f_out.write(line.format(rho_w,mu_w))
                else:
                    f_out.write(line)
                i_line += 1
            



