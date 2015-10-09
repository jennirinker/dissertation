"""
Python dictionary for WindPACT 0.75 MW model from 0.75A08V00 model found in
Archive files
"""
import numpy as np
import scipy.io as scio

fcsvnames = ['S818_2702_geometry.csv',
             'S825_2102_geometry.csv',
             'S826_1602_geometry.csv']
fmatname  = 'WindPACT_geometries'

outdict = {}
for i_file in range(len(fcsvnames)):
    fcsvname = fcsvnames[i_file]
    key      = fcsvname[:9]
    
    BlGeo = np.empty((58,2))
    with open(fcsvname,'rb') as f:
        i_line = 0
        for line in f:
            if i_line > 0:
                BlGeo[i_line-1,:] = [float(x) for x in line.split(',')]
            i_line += 1
          
    outdict[key] = BlGeo
            
        
scio.savemat(fmatname,outdict)


