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

# wind datasets
datasets = ['NREL','fluela','PM06','texastech']
#datasets = ['fluela']

# define turbine name and run name
TurbName = 'WP5.0A04V00'

zRef  = 90.
shear = 0.2
Uref_lo,Uref_hi = 5,22

nbins = 50
NumAvg = 20

cov = 0.1

NumSamps = 20*365*24*6              # number of 10-minute samples

parameters = [['DEL-h','RootMFlp1','MN-m',1000.,12],
              ['DEL-h','HSShftTq','kN-m',1,5],
              ['DEL-h','TwrBsMyt','MN-m',1000,5]]
WindParmSamp = ['Mean_Wind_Speed','Sigma_u','Tau_u','Concentration_u']

DmgDictDir  = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                'processed_data'

# base directory where the stats are stored
BaseDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                'processed_data'
BaseStatDir = os.path.join(BaseDir,'proc_stats')
SaveDir = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
            '2016-02-15_dissertation\\figures'
BaseTurbDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
                'dissertation\\FAST_models\\FAST7'
RSMDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
                'dissertation\\fast_analysis\\fitting_metamodels\\RSMs'

# -----------------------------------------------------------------------------

# load RSM dictionary for that turbine
TurbRSMDictName = '{:s}_RSM.bdat'.format(TurbName)
TurbRSMDictPath = os.path.join(RSMDir,TurbRSMDictName)
with open(TurbRSMDictPath,'rb') as f:
    TurbRSMDict = pickle.load(f)
    
# load turbine dictionary to get hub height
TurbDictPath = os.path.join(BaseTurbDir,TurbName,'parameters',
                      '{:s}_Dict.dat'.format(TurbName))
with open(TurbDictPath,'r') as DictFile:
    zHub = json.load(DictFile)['HH']
    
for dataset in datasets:
    print('Processing dataset \"{:s}\"'.format(dataset))
    
    DmgDictName = 'UltIncr_{:s}.txt'.format(dataset)
    
    # load composite distribution information
    dist_fname = '{:s}_6dist_comp_parms.txt'.format(dataset)
    dist_fpath = os.path.join(BaseDir,dist_fname)
    with open(dist_fpath,'r') as f:
        dist_dict  = json.load(f)
    p_parms_opt = dist_dict['p_parms_opt']
    parms       = dist_dict['parms']
    parms[2] = 'Tau_u'
    
    # get height index closest to hub height
    heights = jr.datasetSpecs(dataset)['IDs']
    iH      = np.abs(heights - zHub).argmin()
    zSamp   = heights[iH]
    if dataset == 'PM06': zSamp = 1.5           # measurement height for Plaine Morte
    
    # calculate wind speed range, sample from empirical distribution
    Usamp_lo,Usamp_hi = Uref_lo*(zSamp/zRef)**shear, \
                    Uref_hi*(zSamp/zRef)**shear
    WindParms = jr.SampleWindParameters(NumSamps,dataset,BaseDir,WindParmSamp,iH,
                             URange=[Usamp_lo,Usamp_hi])
    
    # create wind parameter array
    x = np.empty((NumSamps,4))
    x[:,0],x[:,1],x[:,2],x[:,3] = WindParms[:,0]/(zSamp/zRef)**shear,\
                                  WindParms[:,1]/(WindParms[:,0]/(zSamp/zRef)**shear),\
                                  np.log10(WindParms[:,0]*WindParms[:,2]),\
                                  WindParms[:,3]
    
    del WindParms
    
    # loop through statistics
    iPlot = 0
    
    Fs = np.empty((len(parameters),2))
    ns = np.empty((len(parameters),2,nbins))
    bs = np.empty((len(parameters),2,nbins+1))
    for istat in range(len(parameters)):
        stat,parm,units,scale,m = parameters[istat]
        
        print('    {:s} {:s}'.format(stat,parm))
        
        # load the RSM data for that statistic
        DictKey = '{:s}_{:s}'.format(parm,stat)
        RSMDict = TurbRSMDict[DictKey]
        ps, cs  = RSMDict['ps_red'], RSMDict['cs_red']
            
        # -------------- with temporal coherence ----------------------------------
            
        # calculate mean loads
        Xv = jr.myvander(x,ps)
        mean_loads = np.dot(Xv,cs)
        
        # add randomness
        randn = np.random.normal(size=NumSamps)
        loads = mean_loads + mean_loads*cov*randn
        
        # calculate and save lifetime metric
        Fs[istat,1] = np.mean(np.sort(loads)[-NumAvg:])
        
        # calculate and save histogram
        n, b = np.histogram(loads,bins=nbins,normed=True)
        ns[istat,1] = n
        bs[istat,1] = b
            
        # -------------- without temporal coherence -------------------------------
            
        x[:,3] = 0.
            
        # calculate mean loads
        Xv = jr.myvander(x,ps)
        mean_loads = np.dot(Xv,cs)
        
        # add randomness
        randn = np.random.normal(size=NumSamps)
        loads = mean_loads + mean_loads*cov*randn
        
        # calculate and save lifetime metric
        Fs[istat,0] = np.mean(np.sort(loads)[-NumAvg:])
        
        # calculate and save histogram
        n, b = np.histogram(loads,bins=nbins,normed=True)
        ns[istat,0] = n
        bs[istat,0] = b
        
        del(randn)
                    
    # save data dictionary using pickle
    OutDict = {}
    OutDict['TurbName']   = TurbName
    OutDict['parameters'] = parameters
    OutDict['dataset']    = dataset
    OutDict['Fs']         = Fs
    OutDict['ns']         = ns
    OutDict['bs']         = bs
    DmgDictPath = os.path.join(DmgDictDir,DmgDictName)
    with open(DmgDictPath,'wb') as DictFile:
        pickle.dump(OutDict,DictFile)
    print('\nDictionary {:s} saved.'.format(DmgDictName))


