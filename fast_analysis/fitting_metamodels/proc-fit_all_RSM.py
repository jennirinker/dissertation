"""
Fit response surfaces to all FAST outputs of interest for all wind turbines
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
import os, pickle
import numpy as np
import statsmodels.api as sm


# define turbine name and run name
TurbNames = ['WP0.75A08V00','WP1.5A08V03',
              'WP3.0A02V02','WP5.0A04V00']
#TurbNames = ['WP0.75A08V00']
#RunNames = ['Peregrine','TestRun','SmallRun','BigRun2']
RunName  = 'BigRun2'

FigNum = 1

# statistic and value to fit RSM to
parameters = [['max','RootMFlp1','MN-m',1000.],
              ['DEL-h','RootMFlp1','MN-m',1000.],
              ['max','HSShftTq','kN-m',1],
              ['DEL-h','HSShftTq','kN-m',1],
              ['max','TwrBsMyt','MN-m',1000],
              ['DEL-h','TwrBsMyt','MN-m',1000],
              ['mean','GenPwr','MW',1000.]]
#parameters = [['max','RootMFlp1','MN-m',1000.]]

# base directory where the stats are stored
BaseStatDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                'processed_data\\proc_stats'

alpha = 0.05

# -----------------------------------------------------------------------------
def StatsErr(x,y,p_i,
             alpha=0.05):
    """ Sum-of-squared-error for data x and y with polynomial orders p_i
    """
    
    # perform OLS for all data
    ps_all = jr.GetAllPowers(p_i)
    Xv_all = jr.myvander(x,ps_all)
    results = sm.OLS(y, Xv_all).fit()
    
    # extract significant coefficients
    ps_red = ps_all[results.pvalues <= alpha]
    Xv_red = jr.myvander(x,ps_red)
    
    # fit OLS with reduced coefficients
    cs_red = sm.OLS(y, Xv_red).fit().params
        
    # get reduced model
    yhat   = np.dot(Xv_red,cs_red)
    
    # get residuals and sum
    e    = y - yhat
    err  = np.mean(e ** 2)
    
    return err

# loop through turbines
for TurbName in TurbNames:
    
    TurbRSMDict = {}
    TurbRSMDictName = '{:s}_RSM.bdat'.format(TurbName)
    
    print('Wind turbine {:s}'.format(TurbName))
    
    for stat,parm,units,scale in parameters:
        
        TurbRSMDictKey = '{:s}_{:s}'.format(parm,stat)
        RSMDict         = {}

        # get wind parameters for that run
        WindParms = jr.RunName2WindParms(RunName)
        URefs, Is, Ls, rhos, n_dups = WindParms['URefs'],WindParms['Is'], \
                                      WindParms['Ls'],WindParms['rhos'], \
                                      WindParms['n_dups']
        WindParmsList = [URefs,Is,np.log10(Ls),rhos]
        
        # load the stats data
        x, y = jr.LoadFASTStats(RunName,TurbName,stat,parm)
        
        # ================= optimize polynomial coefficients ==================
        
        # parameterize function for data
        ErrFunc = lambda p: StatsErr(x,y,p)
        
        p0 = np.zeros(x.shape[1])
        results = jr.DiscreteOpt(ErrFunc,p0,
                                 verbose=0)
        p_i = results['p_out']
        
        # ============== plot data and polynomial surface =====================
        
        # get significant powers and coefficients
        ps_all   = jr.GetAllPowers(p_i)
        Xv_all   = jr.myvander(x,ps_all)
        results  = sm.OLS(y, Xv_all).fit()
        cs_all   = results.params
        ps_red   = ps_all[results.pvalues <= alpha]
        Xv_red   = jr.myvander(x,ps_red)
        cs_red   = sm.OLS(y, Xv_red).fit().params
        
        yhat     = np.dot(Xv_red,cs_red)
        es_red   = y - yhat
        perr_red = (yhat - y) / y * 100.
        
        print('{:s} {:s} MSE: {:g}'.format(stat,parm,np.mean(es_red**2)))
        
        # ============== save dictionary =====================
        
        RSMDict['pmax_i'] = p_i
        RSMDict['ps']     = ps_red
        RSMDict['cs']     = cs_red
        RSMDict['es']     = es_red
        RSMDict['p_errs'] = perr_red
        
        TurbRSMDict[TurbRSMDictKey] = RSMDict
        
    # save turbine dictionary
    TurbRSMDictPath = os.path.join('RSMs',TurbRSMDictName)
    with open(TurbRSMDictPath,'wb') as DictFile:
        pickle.dump(TurbRSMDict,DictFile,
                    protocol=pickle.HIGHEST_PROTOCOL)
    


