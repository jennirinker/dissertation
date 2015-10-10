"""
testing function to load turbine dictionary and interpolate blade and 
aerodyn parameters
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import json
#import numpy as np



# define turbine name
TName = '0.75A08V00'

# load turbine model
fTDictName = TName + '_Dict.txt'
with open(fTDictName,'r') as f:
    TurbDict = json.load(f)

BldInterp, ADInterp = jr.InterpolateRotorParams(TurbDict)

for i in range(len(BldInterp)):
    print('{:8.4f}{:8.3f}{:8.2f}{:8.2f}{:12.4g}{:12.4g}{:12.4g}{:12.4g}'\
            .format(*BldInterp[i,:]))
    
print('\n')    
    
for i in range(len(ADInterp)):
    print('{:8.4f}{:8.2f}{:12.5f}{:8.3f}{:8.1f}'\
            .format(*ADInterp[i,:]))