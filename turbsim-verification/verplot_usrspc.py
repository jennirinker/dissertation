""" Script to load a turbsim file and create a plot
	 with subplots for verification. """

import pyts.io.main as io
import matplotlib.pyplot as plt
import numpy as np
import csv

import sys
##libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
libpath = '/home/jrinker/git/dissertation/'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
import os

nticks = 4

# name of file to load
##dname   = '3-PDDs/TS/'
#dname,fignum   = '5-specify_allvals/TS/',1
dname,fignum   = '6-allrho_newsamp/TS/',2
##fname,fignum   = '5pts_NoSc',1
fname   = '5pts_Usr'
#dname = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#            'dissertation\\FAST_models\\WindPACT\\code\\linux_bts'
#fname,fignum   = 'WP1.5A08V03_43333',1

# save image in directory?
saveimg = 0

# construct total file path
inpname = os.path.join(dname,fname + '.inp')
outname = os.path.join(dname,fname + '.bts')
spcname = os.path.join(dname,fname + '.spc')

# read files
tsout = io.readModel(outname)
tsin  = jr.readInput_v2(inpname)
with open(spcname,'r') as f_obj:
    stop, i_line, i_f, NumF = 0, 0, 0, 1e6
    while not stop:
        line = f_obj.readline()
        if i_line == 1:
            contents = line.split(';')
            contents[-1] = contents[-1].split()[0]
            parms = [float(x.split('=')[-1]) for x in contents]
            URef,sig_u,sig_v,sig_w,L_u,L_v,L_w,\
                rho_u,rho_v,rho_w,mu_u,mu_v,mu_w = parms
        elif i_line == 3:
            contents = line.split()
            NumF = int(contents[0])
            S_theo = np.empty((NumF,4))
        elif i_line == 4:
            contents = line.split()
            Scale1 = float(contents[0])
        elif i_line == 5:
            contents = line.split()
            Scale2 = float(contents[0])
        elif i_line == 6:
            contents = line.split()
            Scale3 = float(contents[0])
        elif i_f >= NumF:
            stop = 1
        elif i_line > 10:
            S_unsc = np.asarray([float(x) for \
                                        x in line.split()])
            S_theo[i_f,:] = [S_unsc[0],S_unsc[1]*Scale1,
                             S_unsc[2]*Scale2,S_unsc[3]*Scale3]
            i_f += 1
        i_line += 1

# useful values
y = tsout.grid.y;               # y-grid vector
z = tsout.grid.z[::-1];         # bts
#z = tsout.grid.z;         # wnd
[Y, Z] = np.meshgrid(y, z);     # grid arrays
n_t = tsout.uhub.size;          # (grid.n_t is not correct)
n_f = jr.uniqueComponents(n_t);     # unique components counting DC
rsep = min(tsout.grid.dy,\
           tsout.grid.dz);      # coherence sep distance
ZRef = tsin.zref                # refernce height

# ============== TurbSim ==============

# calculate velocity and turb intesn profiles
U  = jr.TurbSimVelProfile(outname)
Ti = tsout.Ti
sig = U*Ti
		
# calculate spatial coherence
f, Coh = jr.TurbSimSpatCoh(outname,rsep)

# calculate hub-height spectra
Suk, Svk, Swk = jr.TurbSimHHPSDs(outname)

# phase differences
dthetau, dthetav, dthetaw = \
    jr.TurbSimHHPDDs(outname)

# ============== Theory ==============
zhub = tsout.grid.zhub
df = 1./(n_t*tsout.dt)
U_theo  = jr.IEC_VelProfile(z,ZRef,URef)
sig_theo = sig_u*np.ones(z.shape)
Coh_theo = jr.IEC_SpatialCoherence(zhub,URef,rsep,f)
f_theo = S_theo[:,0]
Su_theo = S_theo[:,1]
Sv_theo = S_theo[:,2]
Sw_theo = S_theo[:,3]
thetaPlot = np.linspace(-np.pi,np.pi,100)
fWC_u = jr.wrappedCauchyPDF(thetaPlot,rho_u,mu_u)
fWC_v = jr.wrappedCauchyPDF(thetaPlot,rho_v,mu_v)
fWC_w = jr.wrappedCauchyPDF(thetaPlot,rho_w,mu_w)

# print standard deviations
sigHH_theo = sig_theo[0]
sigHH_TS  = np.std(tsout.uhub)
perc_diff = (sigHH_TS-sigHH_theo)/sigHH_theo*100
print('IEC hub-height sigma:     {:.3f}'.format(sigHH_theo))
print('TurbSim hub-height sigma: {:.3f}'.format(sigHH_TS))
print('Percent difference:       {:.1f}%'.format(perc_diff))
print(np.std(tsout.vhub)/sig_v)

# ============== Plot ==============
axwidth = 0.23;
axheight = 0.23;
xedge = [0.08, 0.41, 0.74];
yedge = [0.09, 0.405, 0.72];

fig = plt.figure(fignum,figsize=(9.05,7.14));
fig.clf();
plt.ion()

# velocity profile
ax1 = plt.axes([xedge[0], yedge[-1], axwidth, axheight])
ax1.scatter( U, Z )
ax1.plot( U_theo, z, 'r')
plt.xlabel('Velocity (m/s)')
plt.ylabel('Height (m)')
plt.locator_params(axis = 'x', nbins = nticks)

# turbulence intensity profile
ax4 = plt.axes([xedge[0], yedge[1], axwidth, axheight])
ax4.scatter( sig, Z )
ax4.plot( sig_theo, z, 'r')
plt.xlabel('Turb. Std. Dev. (-)')
##ax4.scatter( Ti, Z )
##ax4.plot( Ti_theo, z, 'r')
##plt.xlabel('Turb. Intens (-)')
plt.ylabel('Height (m)')
plt.ylabel('Height (m)')

# spatial coherence
ax7 = plt.axes([xedge[0], yedge[0], axwidth, axheight]);
ax7.semilogx(f,Coh);
ax7.semilogx(f,Coh_theo,'r')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Spatial Coherence (-)')

# u-spectrum
ax2 = plt.axes([xedge[1], yedge[-1], axwidth, axheight]);
ax2.loglog(f[1:],Suk[1:-1]/df)
ax2.loglog(f_theo,Su_theo,'r')
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD (m^2/s^2/Hz)')
plt.xlim([1e-3,1e1])
plt.title(fname)

# v-spectrum
ax5 =  plt.axes([xedge[1], yedge[1], axwidth, axheight]);
ax5.loglog(f[1:],Svk[1:-1]/df)
ax5.loglog(f_theo,Sv_theo,'r')
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD (m^2/s^2/Hz)')
plt.xlim([1e-3,1e1])

# w-spectrum
ax8 =  plt.axes([xedge[1], yedge[0], axwidth, axheight]);
ax8.loglog(f[1:],Swk[1:-1]/df)
ax8.loglog(f_theo,Sw_theo,'r')
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD (m^2/s^2/Hz)')
plt.xlim([1e-3,1e1])

# u-dtheta
nbins = 20
ax3 =  plt.axes([xedge[2], yedge[-1], axwidth, axheight]);
ax3.hist(dthetau,bins=nbins,normed=True)
ax3.plot(thetaPlot, fWC_u, 'r')
plt.xlabel('Phase Difference/pi')
plt.ylabel('uhub Phase Differences')
plt.xlim([-np.pi,np.pi])

# v-dtheta
ax6 =  plt.axes([xedge[2], yedge[1], axwidth, axheight]);
ax6.hist(dthetav,bins=nbins,normed=True)
ax6.plot(thetaPlot, fWC_v, 'r')
plt.xlabel('Phase Difference/pi')
plt.ylabel('Uhub Phase Differences')
plt.xlim([-np.pi,np.pi])

# w-dtheta
ax9 =  plt.axes([xedge[2], yedge[0], axwidth, axheight]);
ax9.hist(dthetaw,bins=nbins,normed=True)
ax9.plot(thetaPlot, fWC_w, 'r')
plt.xlabel('Phase Difference/pi')
plt.ylabel('whub Phase Differences')
plt.xlim([-np.pi,np.pi])

fig.show()

if saveimg:
    plt.savefig(dname + fname + '.png')
    print fname + ' saved as image'

