"""
Python dictionary for WindPACT 1.5 MW model from 1.5A08V03 model found in
Archive files
"""
import csv
import numpy as np
import scipy.io as scio

fcsvname = '1.5A08V03_mod.csv'
NBlNodes = 21
NADNodes = 15

#RotorRad = 35
#HubRad   = 1.75
#
#airfoils = ['cylinder.dat','s818_2702.dat','s825_2102.dat','s826_2102.dat']
#cols = ['BlFract','AeroCent','StrcTwst','BMassDen','FlpStff','EdgStff',
#        'GJStff','EAStff','NFoil','Chord']
#
#BladeData = np.empty((NBlNodes,len(cols)))
#
#with open(fcsvname,'rb') as f:
#    reader = csv.DictReader(f)
#    i_line = 0
#    for line in reader:
#        for i_col in range(len(cols)):
#            BladeData[i_line,i_col] = line[cols[i_col]]
#        i_line += 1
#        
#fmatname = fcsvname[:-4] + '.mat'
#outdict = {}
#outdict['airfoils'] = airfoils
#outdict['cols']     = cols
#outdict['BladeData'] = BladeData
#outdict['HubRad']    = HubRad
#outdict['RotorRad'] = RotorRad
#scio.savemat(fmatname,outdict)

fmatname = fcsvname[:-4] + '.mat'
struc = scio.loadmat(fmatname,squeeze_me=True)
Airfoils = [s.strip() for s in struc['airfoils']]
cols     = [s.strip() for s in struc['cols']]
BladeData = struc['BladeData']
RotorRad  = struc['RotorRad']
HubRad    = struc['HubRad']


BladeLen = RotorRad - HubRad
DR       = BladeLen/float(NADNodes)

NodeEdges  = np.linspace(0,1,NADNodes+1) 
NodeCntrs = 0.5*(NodeEdges[1:] + NodeEdges[:-1])
RNodes    = NodeCntrs*BladeLen + HubRad
TwstCntrs  = np.interp(NodeCntrs,BladeData[:,cols.index('BlFract')],
                                BladeData[:,cols.index('StrcTwst')])
ChordCntrs = np.interp(NodeCntrs,BladeData[:,cols.index('BlFract')],
                                BladeData[:,cols.index('Chord')])
AFoilCntrs = np.round(np.interp(NodeCntrs,BladeData[:,cols.index('BlFract')],
                                BladeData[:,cols.index('NFoil')]))
         
for i in range(BladeData.shape[0]):
    print('{:7.5f}{:7.3f}{:10.2f}{:12.2f}{:13.4E}{:13.4E}{:13.4E}{:13.4E}'.format(
            *BladeData[i,:]))
         
for i in range(len(Airfoils)):
    print('\"AeroData/{}\"'.format(Airfoils[i]))
                 
for i in range(NADNodes):
    print('{:8.5f}{:7.2f}{:12.5f}{:7.3f}{:3.0f}{:>13s}'.format(
          RNodes[i],TwstCntrs[i],DR,ChordCntrs[i],AFoilCntrs[i],'NOPRINT'))