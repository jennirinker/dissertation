#! /usr/bin/python
"""
Create user spectrum file
"""
import sys
##libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
libpath = '/home/jrinker/git/dissertation/'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np

# define wind parameters
Uhub  = 10.
sig_u = 2.
L_u   = 340.2
rho_u = 0.2
rho_v = rho_u
rho_w = rho_u

# turbine/simulation parameters
Zhub   = 90.
T      = 600.
dt     = 0.05
fname  = '5pts_Usr.spc'

# create unscaled spectral values
fo, ff, df = 0, 20, 0.0005
fs  = np.arange(fo,ff,df)
sig_v = 0.8*sig_u
sig_w = 0.5*sig_u
L_v   = L_u/8.1*2.7
L_w   = L_u/8.1*0.66
Su    = jr.KaimalSpectrum(fs,L_u/Uhub,sig_u)
Sv    = jr.KaimalSpectrum(fs,L_v/Uhub,sig_v)
Sw    = jr.KaimalSpectrum(fs,L_w/Uhub,sig_w)
NumF  = fs.size

# frequenices/spectral values for scaling
df_s, n_t = 1./T, T/dt
fs_s = np.arange(jr.uniqueComponents(\
        n_t))*df_s
Su_s  = jr.KaimalSpectrum(fs_s,L_u/Uhub,sig_u)
Sv_s  = jr.KaimalSpectrum(fs_s,L_v/Uhub,sig_v)
Sw_s  = jr.KaimalSpectrum(fs_s,L_w/Uhub,sig_w)
Suk_s, Svk_s, Swk_s = Su_s*df_s, Sv_s*df_s, Sw_s*df_s
alpha1 = jr.spectralScale(Suk_s,sig_u,n_t)**2
alpha2 = jr.spectralScale(Svk_s,sig_v,n_t)**2
alpha3 = jr.spectralScale(Swk_s,sig_w,n_t)**2
##Scale1,Scale2,Scale3 = 1.,1.,1.
Scale1,Scale2,Scale3 = alpha1,alpha2,alpha3

with open(fname,'w') as f:
    f.write('-------- User-Defined Spectra (Used only with USRINP ' + \
            'spectral model) ------------------------------------\n')
    f.write('- Kaimal: U={:4.1f};s_u={:4.2f};'.format(Uhub,sig_u) + \
            's_v={:4.2f};s_w={:4.2f};L_u={:5.1f};'.format(sig_v,sig_w,L_u) + \
            'L_v={:5.1f};L_w={:5.1f};r_u={:4.2f};'.format(L_v,L_w,rho_u) + \
            'r_v={:4.2f};r_w={:4.2f}    -\n'.format(rho_v,rho_w))
    f.write('------------------------------------------------------' + \
            '---------------------------------------------------\n')
    f.write('{:<8.0f}        NumUSRf        - Number of '.format(NumF) + \
            'Frequencies [dictates how many lines to read from this file]\n')
    f.write('{:<8.3f}        SpecScale1     - scaling '.format(Scale1) + \
            'factor for the input u-component spectrum\n')
    f.write('{:<8.3f}        SpecScale2     - scaling '.format(Scale2) + \
            'factor for the input v-component spectrum\n')
    f.write('{:<8.3f}        SpecScale3     - scaling '.format(Scale3) + \
            'factor for the input w-component spectrum\n')
    f.write('.........................................................' + \
            '................................................\n')
    f.write('Frequency    u-component PSD   v-component PSD ' + \
            '     w-component PSD\n')
    f.write(' (Hz)           (m^2/s)           (m^2/s)             (m^2/s)\n')
    f.write('-----------------------------------------------' + \
            '----------------------------------------------------------\n')
    for i_f in range(NumF):
        f.write('{:>6.4f}        {:>10.6f}        '.format(fs[i_f],Su[i_f]) + \
                '{:>10.6f}        {:>8.6f}\n'.format(Sv[i_f],Sw[i_f]))
##        f.write('{:>6.4f}        {:>10.6f}        '.format(fs[i_f],alpha1*Su[i_f]) + \
##                '{:>10.6f}        {:>8.6f}\n'.format(alpha2*Sv[i_f],alpha3*Sw[i_f]))
    
    
