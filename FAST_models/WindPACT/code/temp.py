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
#TName = 'WP0.75A08V00'
TName = 'WP1.5A08V03'
#TName = 'WP3.0A02V02'
#TName = 'WP5.0A04V00'

# load turbine model
fTDictName = TName + '_Dict.txt'
with open(fTDictName,'r') as f:
    TurbDict = json.load(f)

TowerInterp = jr.InterpolateTowerParams(TurbDict)

  
    
for i in range(len(TowerInterp)):
    print('{:6.4f}{:12.2f}{:14.5E}{:14.3E}{:14.3E}{:14.3E}'\
            .format(*TowerInterp[i,:]))