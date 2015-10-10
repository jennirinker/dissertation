"""
Load turbine model dictionary and create input files
"""
import json
import numpy as np

ModesInp = 1

# define turbine name
TName = '0.75A08V00'

# load turbine model
fTDictName = TName + '_Dict.txt'
with open(fTDictName,'r') as f:
    TurbDict = json.load(f)

if ModesInp:
    RotorRad = TurbDict['Rotor']['RotDiam']/2.
    SSAngVel = TurbDict['Nacelle']['RatedTipSpeed']/2/np.pi/RotorRad*60.
    HubRad   = TurbDict['Rotor']['HubDiam']/2
    