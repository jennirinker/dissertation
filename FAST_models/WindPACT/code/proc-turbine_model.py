"""
Create a python dictionary describing a turbine from a series of text files
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
import json

# define turbine name
#TName = 'WP0.75A08V00'
TName = 'WP1.5A08V03'
#TName = 'WP3.0A02V02'
#TName = 'WP5.0A04V00'

# **** TURBINE PARAMETERS NOT STORED IN TEXT FILES ****
BldStats = np.concatenate(([0.05],[0.07],
                           np.arange(.1,1.01,.05)))
BldEdges = ((BldStats - BldStats[0])/(BldStats[-1] - BldStats[0])).tolist()
ADEdges  = np.linspace(0,1,16).tolist()
TowerEdges = np.linspace(0,1,10).tolist()

# ===== Create turbine dictionary =====
TurbDict = {}

# ===== Create and add rotor dictionary =====
RotorDict = {}

# load rotor data from text file
n_skip = 5                              # specific for text file
fRotName = 'textfiles//'+TName+'_Blade.txt'
with open(fRotName,'r') as f:
    RotDiam        = [float(x) for x in f.readline().strip('\n').split()][0]
    HubDiam        = [float(x) for x in f.readline().strip('\n').split()][0]
    Airfoils       = [s.split('.')[0] for s in f.readline().strip('\n').split()]
    Damping        = [float(x) for x in f.readline().strip('\n').split()]
    BldSchedFields = f.readline().strip('\n').split('\t')
BldSched = np.genfromtxt(fRotName,skip_header=n_skip).tolist()

# load mode shape data from Modes output
fModesName = 'textfiles//'+TName+'_BldModes.mod'
RotModes = np.empty((3,5))
with open(fModesName,'r') as f:
    i_line = 0
    for line in f:
        if ((i_line >= 20) and (i_line <= 24)):
            RotModes[0,i_line-20] = float(line.split()[1])
            RotModes[1,i_line-20] = float(line.split()[2])
        elif ((i_line >= 33) and (i_line <= 37)):
            RotModes[2,i_line-33] = float(line.split()[1])
        i_line += 1
RotModes = RotModes.tolist()


# save values in blade dictionary
RotorDict['Airfoils']       = Airfoils
RotorDict['BldSched']       = BldSched
RotorDict['BldSchedFields'] = BldSchedFields
RotorDict['Damping']        = Damping
RotorDict['RotDiam']        = RotDiam
RotorDict['HubDiam']        = HubDiam
RotorDict['RotModes']       = RotModes
RotorDict['BldEdges']       = BldEdges
RotorDict['ADEdges']        = ADEdges

# add blade dictionary to turbine dictionary
TurbDict['Rotor'] = RotorDict

# ===== Create and add nacelle dictionary =====
NacDict = {}

# load nacelle data from text file
fRotName = 'textfiles//'+TName+'_Nacelle.txt'
with open(fRotName,'r') as f:
    RatedTipSpeed = [float(x) for x in f.readline().strip('\n').split()][0]
    HubHeight     = [float(x) for x in f.readline().strip('\n').split()][0]

# save values in blade dictionary
NacDict['RatedTipSpeed'] = RatedTipSpeed
NacDict['HubHeight']    = HubHeight

# add blade dictionary to turbine dictionary
TurbDict['Nacelle'] = NacDict

# ===== Create and add tower dictionary =====
TowerDict = {}

# load nacelle data from text file
fRotName = 'textfiles//'+TName+'_Tower.txt'
with open(fRotName,'r') as f:
    HHtoTop   = [float(x) for x in f.readline().strip('\n').split()][0]
    TopDiam   = [float(x) for x in f.readline().strip('\n').split()][0]
    TopThick  = [float(x) for x in f.readline().strip('\n').split()][0]
    BaseDiam  = [float(x) for x in f.readline().strip('\n').split()][0]
    BaseThick = [float(x) for x in f.readline().strip('\n').split()][0]
    ParaMass  = [float(x) for x in f.readline().strip('\n').split()][0]
    TowerDens = [float(x) for x in f.readline().strip('\n').split()][0]
    TowerE    = [float(x) for x in f.readline().strip('\n').split()][0]
    TowerG  = [float(x) for x in f.readline().strip('\n').split()][0]

# save values in blade dictionary
TowerDict['HHtoTop']    = HHtoTop
TowerDict['TopDiam']    = TopDiam
TowerDict['TopThick']   = TopThick
TowerDict['BaseDiam']   = BaseDiam
TowerDict['BaseThick']  = BaseThick
TowerDict['ParaMass']   = ParaMass
TowerDict['TowerDens']  = TowerDens
TowerDict['TowerE']     = TowerE
TowerDict['TowerG']     = TowerG
TowerDict['TowerEdges'] = TowerEdges

# add blade dictionary to turbine dictionary
TurbDict['Tower'] = TowerDict

# ===== Save turbine dictionary =====
fTDictName = TName + '_Dict.txt'
with open(fTDictName,'w') as f:
    json.dump(TurbDict,f)
    
    
    
    