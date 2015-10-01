"""
examining mu of different components
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
import scipy.io as scio
import matplotlib.pyplot as plt
import scipy.stats

fname = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
        'processed_data\\NREL-metadata_allheights.mat'
struc = scio.loadmat(fname)

fields = struc['fields']
data   = struc['values'][:,1:]
hts    = struc['heights']

# strip out values of rho less than 0.06
rho_thresh = 0.0

# get parameter indices
i_u = np.char.strip(fields).tolist().index('mu_u')
i_v = np.char.strip(fields).tolist().index('mu_v')
i_w = np.char.strip(fields).tolist().index('mu_w')
i_rhou = np.char.strip(fields).tolist().index('rho_u')
i_rhov = np.char.strip(fields).tolist().index('rho_v')
i_rhow = np.char.strip(fields).tolist().index('rho_w')

# initialize figure
fig1 = plt.figure(1)
plt.clf()
fig2 = plt.figure(2)
plt.clf()
fig3 = plt.figure(3)
plt.clf()

g_lim = [-4,4]
r_lim = [-np.pi,np.pi]

norm = 0

for iH in range(hts.size):
    
    iu = iH*fields.size + i_u
    iv = iH*fields.size + i_v
    iw = iH*fields.size + i_w
    
    x_u = data[:,iu]
    x_v = data[:,iv]
    x_w = data[:,iw]
    rho_u = data[:,i_rhou]
    rho_v = data[:,i_rhov]
    rho_w = data[:,i_rhow]
    
    i_hirho = np.logical_and(np.logical_and(rho_u >= rho_thresh,
                                            rho_v >= rho_thresh),
                            rho_w >= rho_thresh)
    x_u = x_u[i_hirho]
    x_v = x_v[i_hirho]
    x_w = x_w[i_hirho]
    
    g_u   = np.empty(x_u.shape[0])
    g_v   = np.empty(x_v.shape[0])
    g_w   = np.empty(x_w.shape[0])
    F_emp = (1.+np.arange(rho_w.size))/(1.+rho_w.size)

    g_u[np.argsort(x_u)] = scipy.stats.norm.ppf(F_emp)
    g_v[np.argsort(x_v)] = scipy.stats.norm.ppf(F_emp)
    g_w[np.argsort(x_w)] = scipy.stats.norm.ppf(F_emp)
    
    if norm:
            
        ax1 = fig1.add_subplot(3,2,iH+1)
        p1 = scipy.stats.pearsonr(g_u,g_w)[0]
        ax1.plot(g_lim,g_lim,'k:')
        ax1.scatter(g_u,g_w,s=1)
        ax1.text(2,-4,'r = {:.2f}'.format(p1))

        ax2 = fig2.add_subplot(3,2,iH+1)
        p2 = scipy.stats.pearsonr(g_v,g_w)[0]
        ax2.plot(g_lim,g_lim,'k:')
        ax2.scatter(g_v,g_w,s=1)
        ax2.text(2,-4,'r = {:.2f}'.format(p2))
        
        ax3 = fig3.add_subplot(3,2,iH+1)
        p3 = scipy.stats.pearsonr(g_u,g_v)[0]
        ax3.plot(g_lim,g_lim,'k:')
        ax3.scatter(g_u,g_v,s=1)
        ax3.text(2,-4,'r = {:.2f}'.format(p3))
        
    else:
        
        ax1 = fig1.add_subplot(3,2,iH+1)
        p1 = scipy.stats.pearsonr(g_u,g_w)[0]
        ax1.plot(r_lim,r_lim,'k:')
        ax1.scatter(x_u,x_w,s=1)
        ax1.text(0.45,-0.05,'r = {:.2f}'.format(p1))
        
        ax2 = fig2.add_subplot(3,2,iH+1)
        p2 = scipy.stats.pearsonr(g_v,g_w)[0]
        ax2.plot(r_lim,r_lim,'k:')
        ax2.scatter(x_v,x_w,s=1)
        ax2.text(0.45,-0.05,'r = {:.2f}'.format(p2))
        
        ax3 = fig3.add_subplot(3,2,iH+1)
        p3 = scipy.stats.pearsonr(g_u,g_v)[0]
        ax3.plot(r_lim,r_lim,'k:')
        ax3.scatter(x_u,x_v,s=1)
        ax3.text(0.45,-0.05,'r = {:.2f}'.format(p3))
        
    ax1.set_title('{} m'.format(hts[0,iH]))
    ax2.set_title('{} m'.format(hts[0,iH]))
    ax3.set_title('{} m'.format(hts[0,iH]))
        
        
    
fig1.suptitle('Mu_u and Mu_w (rho_thresh={:.2f})'.format(rho_thresh))
fig2.suptitle('Mu_v and Mu_w (rho_thresh={:.2f})'.format(rho_thresh))
fig3.suptitle('Mu_u and Mu_v (rho_thresh={:.2f})'.format(rho_thresh))

fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()

fig1.savefig('Mu1' + norm*'_norm' + '.png')
fig2.savefig('Mu2' + norm*'_norm' + '.png')
fig3.savefig('Mu3' + norm*'_norm' + '.png')
