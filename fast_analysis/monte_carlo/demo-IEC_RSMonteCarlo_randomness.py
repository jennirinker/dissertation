"""
Demonstrate that we need to add in randomness when sampling response surface
so our results aren't weird
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
TurbNames = ['WP5.0A04V00']
RunName  = 'BigRun2'

FigNum = 1
FigSize = (6.0,4.5)
TickBins = 4

savefig = 0

zRef  = 90.
shear = 0.2
Vref  = 50
Iref  = 0.16
Uref_lo,Uref_hi = 5,22

cov = 0.1

NumSamps = 5000

parameters = [['max','RootMFlp1','MN-m',1000.]]
#parameters = [['max','RootMFlp1','MN-m',1000.],
#              ['DEL-h','RootMFlp1','MN-m',1000.],
#              ['max','HSShftTq','kN-m',1],
#              ['DEL-h','HSShftTq','kN-m',1],
#              ['max','TwrBsMyt','MN-m',1000],
#              ['DEL-h','TwrBsMyt','MN-m',1000],
#              ['mean','GenPwr','MW',1000.]]

# base directory where the stats are stored
BaseStatDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                'processed_data\\proc_stats'
SaveDir = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
            '2016-02-15_dissertation\\figures'
BaseTurbDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
                'dissertation\\FAST_models\\FAST7'
RSMDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
                'dissertation\\fast_analysis\\fitting_metamodels\\RSMs'

# -----------------------------------------------------------------------------


rands = np.random.rand(NumSamps)
randn = np.random.normal(size=NumSamps)

rho = 0.0

for TurbName in TurbNames:
    
    print('Turbine {:s}'.format(TurbName))
    
    # load RSM dictionary
    TurbRSMDictName = '{:s}_RSM.bdat'.format(TurbName)
    TurbRSMDictPath = os.path.join(RSMDir,TurbRSMDictName)
    with open(TurbRSMDictPath,'rb') as f:
        TurbRSMDict = pickle.load(f)
    
    # load look up table and turbine dictionry
    SSPath = os.path.join(BaseTurbDir,TurbName,'steady_state',
                          '{:s}_SS.mat'.format(TurbName))
    SSDict = scio.loadmat(SSPath,squeeze_me=True)
    SSData = SSDict['SS']
    SSFields = [s.rstrip() for s in SSDict['Fields']]
    TurbDictPath = os.path.join(BaseTurbDir,TurbName,'parameters',
                          '{:s}_Dict.dat'.format(TurbName))
    with open(TurbDictPath,'r') as DictFile:
        HubHeight = json.load(DictFile)['HH']
    logL = np.log10(jr.IEC_Lambda1(HubHeight)*8.1)
        
    Uhub_lo = Uref_lo*(HubHeight/zRef)**shear
    Uhub_hi = Uref_hi*(HubHeight/zRef)**shear
    F_lo    = 1 - np.exp(-np.pi*(Uhub_lo/0.4/Vref)**2)
    F_hi    = 1 - np.exp(-np.pi*(Uhub_hi/0.4/Vref)**2)
    
    rands   = (F_hi - F_lo)*rands + F_lo
    Uhubs = 0.4*Vref*np.sqrt(-np.log(1-rands)/np.pi)
    sig_us = Iref*(0.75*Uhubs+5.6)

    Urefs = Uhubs/(HubHeight/zRef)**shear
    Is    = sig_us/Urefs
        
    x = np.empty((NumSamps,4))
    x[:,0],x[:,1],x[:,2],x[:,3] = Urefs,Is,logL,rho
    
    for stat,parm,units,scale in parameters:
        
        print('  {:s} {:s}'.format(stat,parm))
        
        # load the RSM data
        DictKey = '{:s}_{:s}'.format(parm,stat)
        RSMDict = TurbRSMDict[DictKey]
        ps, cs  = RSMDict['ps_red'], RSMDict['cs_red']
        
        # calculate mean loads
        Xv = jr.myvander(x,ps)
        mean_loads = np.dot(Xv,cs)
        
        # add randomness
        loads = mean_loads + mean_loads*cov*randn
        
        plt.figure(FigNum,figsize=FigSize)
        plt.clf()
        
        ax = plt.subplot2grid((3,2), (0,0))
        n,b,p = plt.hist(Urefs,bins=50,histtype='step',normed=True)
        plt.xlabel('Reference Wind Speed (m/s)',fontsize='small')
        ax.set_ylim([0,n.max()*1.25])
        plt.locator_params(nbins=5)
        
        ax = plt.subplot2grid((3,2), (0,1))
        n,b,p = plt.hist(Is,bins=50,histtype='step',normed=True)
        plt.xlabel('Turbulence Intensity',fontsize='small')
        ax.set_ylim([0,n.max()*1.25])
        plt.locator_params(nbins=5)
        
        ax = plt.subplot2grid((3,2), (2,0), colspan=2)
        n,b,p = plt.hist(loads,bins=50,histtype='step',normed=True)
        plt.xlabel('{:s} {:s} ({:s})'.format(stat,parm,units),fontsize='small')
        ax.set_ylim([0,n.max()*1.25])
        ax.set_title('Response Surface Loads with 10% Covariance')
        plt.locator_params(nbins=5)
        xlim = ax.get_xlim()
        
        ax = plt.subplot2grid((3,2), (1,0), colspan=2)
        n,b,p = plt.hist(mean_loads,bins=50,histtype='step',normed=True)
        plt.xlabel('{:s} {:s} ({:s})'.format(stat,parm,units),fontsize='small')
        ax.set_ylim([0,n.max()*1.25])
        ax.set_xlim(xlim)
        ax.set_title('Response Surface Loads')
        plt.locator_params(nbins=5)
        
plt.tight_layout()
            
