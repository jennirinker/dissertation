"""
Write/run TurbSim and FAST files for all four wind turbines with varying rho
but same random seeds
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\nwtc_python_tools'
if (libpath not in sys.path): sys.path.append(libpath)
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import jr_fast
import JR_Library.main as jr
import os, json
import random

# define turbine names
TurbNames = ['WP0.75A08V00','WP1.5A08V03',
              'WP3.0A02V02','WP5.0A04V00']
#TurbNames = ['WP0.75A08V00']

BaseTurbDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\FAST7'
TmplDir     = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\' + \
    'nwtc_python_tools\\templates'
BaseFastDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation\\' + \
    'fast_analysis\\run_fast\\TurbComp'

# concentration parameters to simulate
rhos = [0.0,0.1,0.3]

UHub, IRef, L_u = 15., 0.14, 340.2
zRef = 90.
TSRandLo,TSRandHi = -2147483648,2147483647

# -----------------------------------------------------------------------------

sig_u = UHub * IRef
R1s = [random.randint(TSRandLo, TSRandHi),random.randint(TSRandLo, TSRandHi),
       random.randint(TSRandLo, TSRandHi)]
R2s = [random.randint(TSRandLo, TSRandHi),random.randint(TSRandLo, TSRandHi),
       random.randint(TSRandLo, TSRandHi)]

# change directory to turbine directory to run TurbSim/FAST
os.chdir(BaseFastDir)
    
for TurbName in TurbNames:
    
    print('Turbine {:s}'.format(TurbName))
    
    # define turbine directories
    TurbDir = os.path.join(BaseTurbDir,TurbName)
    ModlDir = TurbDir
    FastDir = os.path.join(BaseFastDir,TurbName)
    
    # get turbine hug-height to back out URef
    DictPath = os.path.join(TurbDir,'parameters',
                            '{:s}_Dict.dat'.format(TurbName))
    with open(DictPath,'r') as DictFile:
        HubHeight = json.load(DictFile)['HH']
    URef = UHub / (HubHeight/zRef) ** 0.2
    
    # loop through concentration parameters
    for iRho in range(len(rhos)):
        
        print('  rho {:d}'.format(iRho))
        
        # define parameters
        FileID = str(iRho)
        rho    = rhos[iRho]
        InpName = '{:s}_{:s}.inp'.format(TurbName,FileID)
        R1, R2  = R1s[iRho], R2s[iRho]
        
        # write turbsim input file
        TSDict = jr.MakeTSDict(TurbName,URef,sig_u,L_u,rho,R1,R2)
        jr.WriteTurbSimInputs(InpName,TSDict,TmplDir,FastDir)
        
        # run TurbSim
        os.system('TurbSim.exe {:s}\\{:s}'.format(TurbName,InpName))
        
        # write FAST input file
        FastName = InpName.rstrip('.inp') 
        BtsName  = FastName + '.bts'
        WindPath = os.path.join(FastDir,BtsName)
        jr_fast.WriteFastADOne(TurbName,WindPath,FastName,ModlDir,FastDir)
        
        # run FAST
        os.system('FAST.exe {:s}\\{:s}'.format(TurbName,FastName+'.fst'))
        
        
        