"""
Calculate the decrease in lifetime for IEC caused by temporal coherence
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
import os, pickle, json
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as scio

# plot style
plt.style.use(jr.stylepath('duke_paper'))

# define turbine name and run name
#TurbNames = ['WP0.75A08V00','WP1.5A08V03',
#              'WP3.0A02V02','WP5.0A04V00']
#RunNames = ['Peregrine','TestRun','SmallRun','BigRun2']
TurbName = 'WP5.0A04V00'
RunName  = 'BigRun2'

savefig = 0

zRef  = 90.
shear = 0.2
Vrefs  = [50,42.5,37.5]
Irefs  = [0.16,0.14,0.12]
Uref_lo,Uref_hi = 5,22
rhos = [0,0.1,0.2,0.3]

nbins = 50

cov = 0.1

NumSamps = 20*365*24*6              # number of 10-minute samples

parameters = [['DEL-h','RootMFlp1','MN-m',1000.,12],
              ['DEL-h','HSShftTq','kN-m',1,5],
              ['DEL-h','TwrBsMyt','MN-m',1000,5]]

# base directory where the stats are stored
BaseStatDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                'processed_data\\proc_stats'
SaveDir = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
            '2016-02-15_dissertation\\figures'
BaseTurbDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
                'dissertation\\FAST_models\\FAST7'
RSMDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
                'dissertation\\fast_analysis\\fitting_metamodels\\RSMs'
DmgDictDir  = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                'processed_data'
DmgDictName = 'DmgIncr_IEC.txt'

# -----------------------------------------------------------------------------

# load RSM dictionary for that turbine
TurbRSMDictName = '{:s}_RSM.bdat'.format(TurbName)
TurbRSMDictPath = os.path.join(RSMDir,TurbRSMDictName)
with open(TurbRSMDictPath,'rb') as f:
    TurbRSMDict = pickle.load(f)

# load turbine dictionary
TurbDictPath = os.path.join(BaseTurbDir,TurbName,'parameters',
                      '{:s}_Dict.dat'.format(TurbName))
with open(TurbDictPath,'r') as DictFile:
    HubHeight = json.load(DictFile)['HH']
logL = np.log10(jr.IEC_Lambda1(HubHeight)*8.1)

# loop through reference wind speeds
Fs    =  np.empty((len(parameters),len(rhos),len(Vrefs),len(Irefs)))
ns    = np.empty((len(parameters),len(rhos),len(Vrefs),len(Irefs),nbins))
bs    = np.empty((len(parameters),len(rhos),len(Vrefs),len(Irefs),nbins+1))
for iVref in range(len(Vrefs)):
    Vref = Vrefs[iVref]
    print('Vref {:.1f}'.format(Vref))
    # loop through reference turbulence intensities
    for iIref in range(len(Irefs)):
        Iref = Irefs[iIref]
        print('  Iref {:.1f}'.format(Iref))
    
        # calculate wind speed limits
        Uhub_lo = Uref_lo*(HubHeight/zRef)**shear
        Uhub_hi = Uref_hi*(HubHeight/zRef)**shear
        F_lo    = 1 - np.exp(-np.pi*(Uhub_lo/0.4/Vref)**2)
        F_hi    = 1 - np.exp(-np.pi*(Uhub_hi/0.4/Vref)**2)    
        
        # loop through statistics
        iPlot = 0
        
        for istat in range(len(parameters)):
            stat,parm,units,scale,m = parameters[istat]
            
            print('    {:s} {:s}'.format(stat,parm))
            
            # load the RSM data for that statistic
            DictKey = '{:s}_{:s}'.format(parm,stat)
            RSMDict = TurbRSMDict[DictKey]
            ps, cs  = RSMDict['ps_red'], RSMDict['cs_red']
                
            # sample random numbers
            rands = np.random.rand(NumSamps)
            
            # scale random numbers, get wind speeds and turbulence std devs
            rands   = (F_hi - F_lo)*rands + F_lo
            Uhubs = 0.4*Vref*np.sqrt(-np.log(1-rands)/np.pi)
            sig_us = Iref*(0.75*Uhubs+5.6)
            Urefs = Uhubs/(HubHeight/zRef)**shear
            Is    = sig_us/Urefs
            
            del(rands,Uhubs,sig_us)
            
            # loop through rho values
            ymax = 0
            for iRho in range(len(rhos)):
                rho = rhos[iRho]
                
                print('      rho {:d}'.format(iRho))
                
                x = np.empty((NumSamps,4))
                x[:,0],x[:,1],x[:,2],x[:,3] = Urefs,Is,logL,rho
                
                # calculate mean loads
                Xv = jr.myvander(x,ps)
                mean_loads = np.dot(Xv,cs)
                
                # add randomness
                randn = np.random.normal(size=NumSamps)
                loads = mean_loads + mean_loads*cov*randn
                
                del(randn)
                
                # calculate and save lifetime metric
                Fs[istat,iRho,iVref,iIref] = np.sum(loads**m)
                
                # calculate and save histogram
                n, b = np.histogram(loads,bins=nbins,normed=True)
                ns[istat,iRho,iVref,iIref] = n
                bs[istat,iRho,iVref,iIref] = b
                
# save data dictionary using pickle
OutDict = {}
OutDict['TurbName']   = TurbName
OutDict['parameters'] = parameters
OutDict['Vrefs']      = Vrefs
OutDict['Irefs']      = Irefs
OutDict['rhos']       = rhos
OutDict['Fs']         = Fs
OutDict['ns']         = ns
OutDict['bs']         = bs
DmgDictPath = os.path.join(DmgDictDir,DmgDictName)
with open(DmgDictPath,'wb') as DictFile:
    pickle.dump(OutDict,DictFile)
print('\nResults saved.')


