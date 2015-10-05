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
from mpl_toolkits.mplot3d import Axes3D

fname = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
        'processed_data\\NREL-metadata_allheights.mat'
struc = scio.loadmat(fname)

fields = struc['fields']
data   = struc['values'][:,1:]
hts    = struc['heights']

# strip out values of rho less than 0.06
rho_thresh = 0.06

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

h_lim = [0,2.5]
r_lim = [-np.pi,np.pi]

norm = 0

for iH in range(hts.size):
    
    iu = iH*fields.size + i_u
    iv = iH*fields.size + i_v
    iw = iH*fields.size + i_w
    
    X_u = data[:,iu]
    X_v = data[:,iv]
    X_w = data[:,iw]
    rho_u = data[:,i_rhou]
    rho_v = data[:,i_rhov]
    rho_w = data[:,i_rhow]
    
    i_hirho = np.logical_and(np.logical_and(rho_u >= rho_thresh,
                                            rho_v >= rho_thresh),
                            rho_w >= rho_thresh)
    X_u = X_u[i_hirho]
    X_v = X_v[i_hirho]
    X_w = X_w[i_hirho]
    
    g_u   = np.empty(X_u.shape[0])
    g_v   = np.empty(X_v.shape[0])
    g_w   = np.empty(X_w.shape[0])
    F_emp = (1.+np.arange(rho_w.size))/(1.+rho_w.size)

    g_u[np.argsort(X_u)] = scipy.stats.norm.ppf(F_emp)
    g_v[np.argsort(X_v)] = scipy.stats.norm.ppf(F_emp)
    g_w[np.argsort(X_w)] = scipy.stats.norm.ppf(F_emp)
    
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
        
#        ax1 = fig1.add_subplot(3,2,iH+1,projection='3d')
#        p1 = scipy.stats.pearsonr(g_u,g_w)[0]
#        theta, phi = X_u, X_w
#        x = np.cos(theta)*np.sin(phi)
#        y = np.sin(theta)*np.sin(phi)
#        z = np.cos(phi)
#        ax1.scatter(x,y,z,s=1)
##        ax1.text(0.45,-0.05,'r = {:.2f}'.format(p1))
#        
#        ax2 = fig2.add_subplot(3,2,iH+1,projection='3d')
#        p2 = scipy.stats.pearsonr(g_v,g_w)[0]
#        theta, phi = X_v, X_w
#        x = np.cos(theta)*np.sin(phi)
#        y = np.sin(theta)*np.sin(phi)
#        z = np.cos(phi)
#        ax2.scatter(x,y,z,s=1)
##        ax2.text(0.45,-0.05,'r = {:.2f}'.format(p2))
#        
#        ax3 = fig3.add_subplot(3,2,iH+1,projection='3d')
#        p3 = scipy.stats.pearsonr(g_u,g_v)[0]
#        theta, phi = X_u, X_v
#        x = np.cos(theta)*np.sin(phi)
#        y = np.sin(theta)*np.sin(phi)
#        z = np.cos(phi)
#        ax3.scatter(x,y,z,s=1)
##        ax3.text(0.45,-0.05,'r = {:.2f}'.format(p3))
        
#        ax1 = fig1.add_subplot(3,2,iH+1)
#        p1 = scipy.stats.pearsonr(g_u,g_w)[0]
#        ax1.plot(r_lim,r_lim,'k:')
#        ax1.scatter(X_u,X_w,s=1)
#        ax1.text(0.45,-0.05,'r = {:.2f}'.format(p1))
#        
#        ax2 = fig2.add_subplot(3,2,iH+1)
#        p2 = scipy.stats.pearsonr(g_v,g_w)[0]
#        ax2.plot(r_lim,r_lim,'k:')
#        ax2.scatter(X_v,X_w,s=1)
#        ax2.text(0.45,-0.05,'r = {:.2f}'.format(p2))
#        
#        ax3 = fig3.add_subplot(3,2,iH+1)
#        p3 = scipy.stats.pearsonr(g_u,g_v)[0]
#        ax3.plot(r_lim,r_lim,'k:')
#        ax3.scatter(X_u,X_v,s=1)
#        ax3.text(0.45,-0.05,'r = {:.2f}'.format(p3))
        
        ax1 = fig1.add_subplot(3,2,iH+1)
        dtheta = jr.wrap(X_u-X_w)
        r      = jr.samplePhaseCoherence(dtheta)[0]
        ax1.hist(dtheta,bins=20,histtype='step',normed=True)
        ax1.text(0.7,0.9,'r = {:.2f}'.format(r),transform=ax1.transAxes)
        ax1.set_ylim(h_lim)
        ax1.set_xlim([-np.pi,np.pi])
        
        ax2 = fig2.add_subplot(3,2,iH+1)
        dtheta = jr.wrap(X_v-X_w)
        r      = jr.samplePhaseCoherence(dtheta)[0]
        ax2.hist(dtheta,bins=20,histtype='step',normed=True)
        ax2.text(0.7,0.9,'r = {:.2f}'.format(r),transform=ax2.transAxes)
        ax2.set_ylim(h_lim)
        ax2.set_xlim([-np.pi,np.pi])
        
        ax3 = fig3.add_subplot(3,2,iH+1)
        dtheta = jr.wrap(X_u-X_v)
        r      = jr.samplePhaseCoherence(dtheta)[0]
        ax3.hist(dtheta,bins=20,histtype='step',normed=True)
        ax3.text(0.7,0.9,'r = {:.2f}'.format(r),transform=ax3.transAxes)
        ax3.set_ylim(h_lim)
        ax3.set_xlim([-np.pi,np.pi])
        
    ax1.set_title('{} m'.format(hts[0,iH]))
    ax2.set_title('{} m'.format(hts[0,iH]))
    ax3.set_title('{} m'.format(hts[0,iH]))
    
    ax1.set_xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
    ax1.set_xticklabels(['-1.0','-0.5','0.0','0.5','1.0'])
    ax2.set_xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
    ax2.set_xticklabels(['-1.0','-0.5','0.0','0.5','1.0'])
    ax3.set_xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
    ax3.set_xticklabels(['-1.0','-0.5','0.0','0.5','1.0'])
        
        
    
fig1.suptitle('Mu_u and Mu_w (rho_thresh={:.2f})'.format(rho_thresh))
fig2.suptitle('Mu_v and Mu_w (rho_thresh={:.2f})'.format(rho_thresh))
fig3.suptitle('Mu_u and Mu_v (rho_thresh={:.2f})'.format(rho_thresh))

fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()

#fig1.savefig('Mu1' + norm*'_norm' + '.png')
#fig2.savefig('Mu2' + norm*'_norm' + '.png')
#fig3.savefig('Mu3' + norm*'_norm' + '.png')
