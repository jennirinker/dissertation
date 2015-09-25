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

Vhub   = 10.
Zhub   = 90.
Tclass = 'B'
T      = 600.
dt     = 0.05
fname  = 'UsrSpec.inp'

##fo, ff, df = 0, 20, 0.001
##freqs  = np.arange(fo,ff,df)
df = 1./T
freqs = np.arange(6001)*df
Su,Sv,Sw = jr.IEC_PSDs(Zhub,Vhub,Tclass,freqs)
NumF   = freqs.size

n_t    = NumF*2 + 1                 # number of time steps
Iref   = jr.IEC_Iref(Tclass)        # reference turbulence intensity
sigma1 = jr.IEC_Sigma1(Iref,Vhub)   # longitudinal standard deviation
sigma2 = 0.8*sigma1                 # lateral standard deviaiton
sigma3 = 0.5*sigma1                 # vertical standard deviation
Suk, Svk, Swk = Su*df, Sv*df, Sw*df
alpha1 = 1./jr.spectralScale(Suk,sigma1,n_t)**2
alpha2 = 1./jr.spectralScale(Svk,sigma2,n_t)**2
alpha3 = 1./jr.spectralScale(Swk,sigma3,n_t)**2
Scale1,Scale2,Scale3 = 1.,1.,1.

with open(fname,'w') as f:
    f.write('-------- User-Defined Spectra (Used only with USRINP ' + \
            'spectral model) ------------------------------------\n')
    f.write('-        The Kaimal spectra IEC 61400-1 Ed. 3 for Vhub' + \
            '={:.0f} m/s; Zhub={:.0f} m; Class="{:s}";'.format(Vhub,Zhub,Tclass) \
            + '                    -\n')
    f.write('------------------------------------------------------' + \
            '---------------------------------------------------\n')
    f.write('{:<8.0f}        NumUSRf        - Number of '.format(NumF) + \
            'Frequencies [dictates how many lines to read from this file]\n')
    f.write('{:<8.2f}        SpecScale1     - scaling '.format(Scale1) + \
            'factor for the input u-component spectrum\n')
    f.write('{:<8.2f}        SpecScale2     - scaling '.format(Scale2) + \
            'factor for the input v-component spectrum\n')
    f.write('{:<8.2f}        SpecScale3     - scaling '.format(Scale3) + \
            'factor for the input w-component spectrum\n')
    f.write('.........................................................' + \
            '................................................\n')
    f.write('Frequency    u-component PSD   v-component PSD ' + \
            '     w-component PSD\n')
    f.write(' (Hz)           (m^2/s)           (m^2/s)             (m^2/s)\n')
    f.write('-----------------------------------------------' + \
            '----------------------------------------------------------\n')
    for i_f in range(NumF):
##        f.write('{:>6.3f}        {:>10.6f}        '.format(freqs[i_f],Su[i_f]) + \
##                '{:>10.6f}        {:>8.6f}\n'.format(Sv[i_f],Sw[i_f]))
        f.write('{:>6.3f}        {:>10.6f}        '.format(freqs[i_f],alpha1*Su[i_f]) + \
                '{:>10.6f}        {:>8.6f}\n'.format(alpha2*Sv[i_f],alpha3*Sw[i_f]))
    
    
