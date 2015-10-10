"""
testing function to load turbine dictionary and interpolate blade and 
aerodyn parameters
"""
import json
#import numpy as np

def InterpolateRotorParams(TurbDict): 
    """ Interpolate blade structural properties and aerodynamic properties
    
        Args:
            TurbDict (dictionary): turbine properties
            
        Returns:
            BldInterp (numpy array): interpolated blade properties
            ADInterp (numpy array): interpolated aerodynamic properties
    """
    import numpy as np
    
    # extract information from turbine dictionary
    BldSched = np.array(TurbDict['Rotor']['BldSched'])
    Fields   = TurbDict['Rotor']['BldSchedFields']
    BldEdges = np.array(TurbDict['Rotor']['BldEdges'])
    ADEdges  = np.array(TurbDict['Rotor']['ADEdges'])
    RotDiam  = TurbDict['Rotor']['RotDiam']
    HubDiam  = TurbDict['Rotor']['HubDiam']

    # place structural values in array
    GenAxLoc = BldSched[:,Fields.index('Gen. Axis Loc.')]
    StrcTwst = BldSched[:,Fields.index('Twist')]
    BMassDen = BldSched[:,Fields.index('Unit Weight')]
    FlpStff  = BldSched[:,Fields.index('EIFlap')]
    EdgStff  = BldSched[:,Fields.index('EIEdge')]
    GJStff   = BldSched[:,Fields.index('GJ')]
    EAStff   = BldSched[:,Fields.index('EA')] 
    Chord    = BldSched[:,Fields.index('Chord')] 
    AFID     = BldSched[:,Fields.index('Airfoil ID')] 
    
    # calculate intermediate parameters
    Station  = np.array(BldSched[:,Fields.index('Station')])
    BlFract  = (Station - Station[0])/(Station[-1] - Station[0])
    ACOff    = Chord*(0.25 - GenAxLoc)
    BladeLen = RotDiam/2. - HubDiam/2.
    ENodes = 0.5*(ADEdges[:-1] + ADEdges[1:])
    RNodes = ENodes*BladeLen + HubDiam/2.
    
    # interpolate blade structural parameters
    ChordInt = np.interp(BldEdges,BlFract,Chord)
    AeroCentInt = np.interp(BldEdges,BlFract,ACOff)/ChordInt + 0.25
    StrcTwstInt = np.interp(BldEdges,BlFract,StrcTwst)
    BMassDenInt = np.interp(BldEdges,BlFract,BMassDen)
    FlpStffInt  = np.interp(BldEdges,BlFract,FlpStff)
    EdgStffInt  = np.interp(BldEdges,BlFract,EdgStff)
    GJStffInt   = np.interp(BldEdges,BlFract,GJStff)
    EAStffInt   = np.interp(BldEdges,BlFract,EAStff)
    
    # place parameters in array
    BldInterp = np.empty((BldEdges.size,8))
    BldInterp[:,0] = BldEdges
    BldInterp[:,1] = AeroCentInt
    BldInterp[:,2] = StrcTwstInt
    BldInterp[:,3] = BMassDenInt
    BldInterp[:,4] = FlpStffInt
    BldInterp[:,5] = EdgStffInt
    BldInterp[:,6] = GJStffInt
    BldInterp[:,7] = EAStffInt
    
    # interpolate aerodynamic parameters
    TwstInt  = np.interp(ENodes,BlFract,StrcTwst)
    ChordInt = np.interp(ENodes,BlFract,Chord)
    AFInt    = np.round(np.interp(ENodes,BlFract,AFID))
    AFInt[1] = 2
    
    # place parameters in array
    ADInterp = np.empty((ADEdges.size-1,5))
    ADInterp[:,0] = RNodes
    ADInterp[:,1] = TwstInt
    ADInterp[:,2] = BladeLen/(ADEdges.size - 1)
    ADInterp[:,3] = ChordInt
    ADInterp[:,4] = AFInt
    print(ENodes)
    print(BlFract)
    print(Chord)
    
    return BldInterp, ADInterp

# define turbine name
TName = '0.75A08V00'

# load turbine model
fTDictName = TName + '_Dict.txt'
with open(fTDictName,'r') as f:
    TurbDict = json.load(f)

BldInterp, ADInterp = InterpolateRotorParams(TurbDict)

for i in range(len(BldInterp)):
    print('{:8.4f}{:8.3f}{:8.2f}{:8.2f}{:12.4g}{:12.4g}{:12.4g}{:12.4g}'\
            .format(*BldInterp[i,:]))
    
print('\n')    
    
for i in range(len(ADInterp)):
    print('{:8.4f}{:8.2f}{:12.5f}{:8.3f}{:8.1f}'\
            .format(*ADInterp[i,:]))