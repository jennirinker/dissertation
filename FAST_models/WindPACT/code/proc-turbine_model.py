"""
Create a python dictionary describing a turbine from a series of text files
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
import json

TName = '0.75A08V00'

# **** TURBINE PARAMETERS NOT STORED IN TEXT FILES ****
BldStats = np.concatenate(([0.05],[0.07],
                           np.arange(.1,1.01,.05)))
BldEdges = ((BldStats - BldStats[0])/(BldStats[-1] - BldStats[0])).tolist()
ADEdges  = np.linspace(0,1,16).tolist()

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
RotModes = np.zeros((3,5)).tolist()

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

# load rotor data from text file
fRotName = 'textfiles//'+TName+'_Nacelle.txt'
with open(fRotName,'r') as f:
    RatedTipSpeed = [float(x) for x in f.readline().strip('\n').split()][0]

# save values in blade dictionary
NacDict['RatedTipSpeed'] = RatedTipSpeed

# add blade dictionary to turbine dictionary
TurbDict['Nacelle'] = NacDict

# ===== Save turbine dictionary =====
fTDictName = TName + '_Dict.txt'
with open(fTDictName,'w') as f:
    json.dump(TurbDict,f)