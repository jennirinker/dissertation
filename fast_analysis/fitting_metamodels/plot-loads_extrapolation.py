"""
Plot the distributions of a load of interest for the fine run
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.stats

# Duke style
plt.style.use(jr.stylepath('duke_paper'))
#plt.style.use(jr.stylepath('duke_presentation'))

# turbine and run name
TurbName = 'WP5.0A04V00'
RunName  = 'Fine'

#FigNum,dist_name = 1, 'exponweib'
#FigNum,dist_name = 2, 'gumbel_r'
FigNum,dist_name = 3, 'genextreme'
savefig = 0
#FigSize = (6.0,8.0)
FigSize = (8.0,6.0) # landscape figure
#FigSize = (6.0,4.5) # vertical half-page

# statistic and value to fit RSM to
#stat,parm,units,scale = 'max','RootMFlp1','MN-m',1000.
#stat,parm,units,scale = 'DEL-h','RootMFlp1','MN-m',1000.
#stat,parm,units,scale = 'max','HSShftTq','kN-m',1
#stat,parm,units,scale = 'DEL-h','HSShftTq','kN-m',1
#stat,parm,units,scale = 'max','TwrBsMyt','MN-m',1000
#stat,parm,units,scale = 'DEL-h','TwrBsMyt','MN-m',1000
parameters = [['max','RootMFlp1','MN-m',1000.],
              ['DEL-h','RootMFlp1','MN-m',1000.],
              ['max','HSShftTq','kN-m',1],
              ['DEL-h','HSShftTq','kN-m',1],
              ['max','TwrBsMyt','MN-m',1000],
              ['DEL-h','TwrBsMyt','MN-m',1000]]
#parameters = [['max','RootMFlp1','MN-m',1000.],
#              ['DEL-h','RootMFlp1','MN-m',1000.],
#              ['max','TwrBsMyt','MN-m',1000],
#              ['DEL-h','TwrBsMyt','MN-m',1000]]

#dist_name = 'genextreme'
#dist_name = 'gumbel_r'
#plot_type = 'pdf'
#plot_type = 'cdf'
plot_type = 'res'
Q = 1

# base directory where the stats are stored
BaseStatDir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                'processed_data\\proc_stats'
SaveDir = 'C:\\Users\\jrinker\\Dropbox\\presentations\\' + \
            '2016-01-18_NWTC\\figures'

colors = ['#235F9C', '#C0504D', '#F79646']

# -----------------------------------------------------------------------------

# get wind parameters for that run
WindParms = jr.RunName2WindParms(RunName)
URefs, Is, Ls, rhos, n_dups = WindParms['URefs'],WindParms['Is'], \
                              WindParms['Ls'],WindParms['rhos'], \
                              WindParms['n_dups']
WindParms = [URefs,Is,np.log10(Ls),rhos]
WindParmStr = ['$U$','$\sigma_u$','log$_{10}(L)$',r'$\rho$']


# ----------------------------------------------------------------------------

# calculate exceednace probability for 50-year return period
t_years = 50.
t_10min = t_years*365*24*6
p_10min = 1./t_10min

print('Results for {:s} with Q = {:.2f}'.format(dist_name,Q))

# initialize figure
fig = plt.figure(FigNum,figsize=FigSize)
plt.clf()
       
fexts = np.empty((len(parameters),len(WindParms[0]),len(WindParms[3])))
for istat in range(len(parameters)):
    
    stat,parm,units,scale = parameters[istat]
    FigTitle = 'plot-fine_load_dists_{:s}_{:s}.png'.format(parm,stat)

    x, y = jr.LoadFASTStats(RunName,TurbName,stat,parm,
                            scale=scale)
    
    titles = ['Below Rated','Above Rated']
    for iU in range(len(WindParms[0])):
        
        # mask data for that wind speed
        U = WindParms[0][iU]
        x_mask = x[x[:,0] == U,:]
        y_mask = y[x[:,0] == U]
        
        # create axes
        ax = plt.subplot(3,4,2*istat+1+iU)
        
        # loop through rhos, fitting and plotting
        ymin,ymax = 1e10,-1e10
        bins = np.linspace(np.floor(y_mask.min()),np.ceil(y_mask.max()),21)
        for iRho in range(len(WindParms[3])):
            
            # mask data for that rho value
            rho    = WindParms[3][iRho]
            y_data = y_mask[x_mask[:,3] == rho]
            x_data = x_mask[x_mask[:,3] == rho]
            F      = (np.arange(y_data.size) +  1.) / (y_data.size + 1.)
            
            # fit distribution
            if Q == 1:
                y_T = float('inf')
            else:
                y_T = np.sort(y_data)[Q*y_data.size]
            p_fit,pGP_fit = jr.fitcompositeparameters(y_data,dist_name,
                                                      x_T=y_T)
                                                      
            
            # calculate residuals, print average residual magnitude
            residuals = jr.compositeCDF(np.sort(y_data),dist_name,p_fit,
                                        x_T=y_T,p_GP=pGP_fit) - F
            mean_res  = np.mean(np.abs(residuals))
            print('Mean residual magnitude (rho{:d}): {:.3f}'.format(iRho,mean_res))
            
            # calculate extrapolated load
            fext = jr.compositeISF(p_10min,dist_name,p_fit,
                                   x_T=y_T,p_GP=pGP_fit)
            fexts[istat,iU,iRho] = fext
            
            if plot_type == 'pdf':
                yPlot  = np.linspace(bins[0],bins[-1]*1.1,200)
#                pdf_fit = dist.pdf(yPlot,*p_fit[:-2],
#                                   loc=p_fit[-2],scale=p_fit[-1])
                pdf_fit = jr.compositePDF(yPlot,dist_name,p_fit,
                                          x_T=y_T,p_GP=pGP_fit)
                                   
                n, bins, patches = plt.hist(y_data,bins=bins,
                     normed=True,histtype='step',label=r'$\rho$ = {:.1f}'.format(rho))
                plt.plot(yPlot,pdf_fit,color=colors[iRho])
                ymin = 0
                ymax = max(ymax,pdf_fit.max()*1.1)
                
                plt.locator_params(nbins=4)
                
            elif plot_type == 'logisf':
                yPlot = np.linspace(y_data.min(),y_data.max(),201)
#                cdf_fit = dist.cdf(yPlot,*p_fit[:-2],
#                                   loc=p_fit[-2],scale=p_fit[-1])
                cdf_fit = jr.compositeCDF(yPlot,dist_name,p_fit,
                                   x_T=y_T,p_GP=pGP_fit)
                
                plt.semilogy(np.sort(y_data),1-F)
                plt.semilogy(yPlot,1-cdf_fit,':',
                             color=colors[iRho],marker='o',mec=colors[iRho],
                             mfc='none',markevery=20)
                ymin = 1e-3
                ymax = 1
                
            elif plot_type == 'cdf':
                yPlot = np.linspace(y_data.min(),y_data.max(),201)
#                cdf_fit = dist.cdf(yPlot,*p_fit[:-2],
#                                   loc=p_fit[-2],scale=p_fit[-1])
                cdf_fit = jr.compositeCDF(yPlot,dist_name,p_fit,
                                   x_T=y_T,p_GP=pGP_fit)
                
                plt.plot(np.sort(y_data),1-F)
                plt.plot(yPlot,1-cdf_fit,':',
                             color=colors[iRho],marker='o',mec=colors[iRho],
                             mfc='none',markevery=20)
                ymin = 0
                ymax = 1
                
            elif plot_type == 'res':
                cdf_fit = jr.compositeCDF(np.sort(y_data),dist_name,p_fit,
                                   x_T=y_T,p_GP=pGP_fit)
                residuals = cdf_fit - F
                
                plt.plot(np.sort(y_data),residuals)
                
                ymin = min(ymin,residuals.min()*1.4)
                ymax = max(ymax,residuals.max()*1.4)
                     
            
        ax.set_ylim([ymin,ymax])
        
        ax.set_title(titles[iU],fontsize='small')
        ax.set_xlabel('{:s} {:s}\n({:s})'.format(stat,parm,units),fontsize='small')
    
#    plt.legend() 
    plt.tight_layout()
    
    if savefig:
        SavePath = os.path.join(SaveDir,FigTitle)
        fig.savefig(SavePath)
        print('Figure saved.')
        
    print('Extrap results {:s} {:s}'.format(stat,parm))
    for iU in range(2):
        extraps = fexts[istat,iU,:]
        perc_incr = 100.*(extraps - extraps[0])/extraps[0]
#        print('{:8.1f} {:8.1f} {:8.1f}'.format(*fexts[istat,iU,:]))
        print('{:8.1f}% {:8.1f}%'.format(*perc_incr[1:]))
        
