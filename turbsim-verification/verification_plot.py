""" Script to load a turbsim file and create a plot
	 with subplots for verification. """

import pyts.io.main as io
import matplotlib.pyplot as plt
import numpy as np

import sys
##libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
libpath = '/home/jrinker/git/dissertation/'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr


# name of file to load
dname   = '3-PDDs/TS/'
fname,fignum   = '5pts_NoSc',1
##fname,fignum   = '5pts_Usr',2

# hard-coded temporal coherence parameters
rho = 0.2           # concentration parameter
mu = 0.        # location parameter

# save image in directory?
saveimg = 0

# construct total file path
inpname = dname + fname + '.inp';
outname = dname + fname + '.wnd';

# read file
tsout = io.readModel(outname);
tsin  = jr.readInput_v2(inpname);

# useful values
y = tsout.grid.y;               # y-grid vector
# z = tsout.grid.z[::-1];         # bts
z = tsout.grid.z;         # wnd
[Y, Z] = np.meshgrid(y, z);     # grid arrays
n_t = tsout.uhub.size;          # (grid.n_t is not correct)
n_f = jr.uniqueComponents(n_t);     # unique components counting DC
rsep = min(tsout.grid.dy,\
           tsout.grid.dz);      # coherence sep distance

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
Vhub = tsout.UHUB
turbc = tsin.turbc
df = 1./(n_t*tsout.dt)
U_IEC  = jr.IEC_VelProfile(z,zhub,Vhub)
Ti_IEC = jr.IEC_TiProfile(z,zhub,Vhub,turbc)
sig_IEC = U_IEC*Ti_IEC
Coh_IEC = jr.IEC_SpatialCoherence(zhub,Vhub,rsep,f);
Su_IEC, Sv_IEC, Sw_IEC = jr.IEC_PSDs(zhub,Vhub,turbc,f);
Suk_IEC = Su_IEC*df
Svk_IEC = Sv_IEC*df
Swk_IEC = Sw_IEC*df
thetaPlot = np.linspace(-np.pi,np.pi,100)
fWC = jr.wrappedCauchyPDF(thetaPlot,rho,mu)

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
ax1.plot( U_IEC, z, 'r')
plt.xlabel('Velocity (m/s)')
plt.ylabel('Height (m)')

# turbulence intensity profile
ax4 = plt.axes([xedge[0], yedge[1], axwidth, axheight])
ax4.scatter( sig, Z )
ax4.plot( sig_IEC, z, 'r')
plt.xlabel('Turb. Intens (-)')
##ax4.scatter( Ti, Z )
##ax4.plot( Ti_IEC, z, 'r')
##plt.xlabel('Turb. Intens (-)')
plt.ylabel('Height (m)')

# spatial coherence
ax7 = plt.axes([xedge[0], yedge[0], axwidth, axheight]);
ax7.semilogx(f,Coh);
ax7.semilogx(f,Coh_IEC,'r')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Spatial Coherence (-)')

# u-spectrum
ax2 = plt.axes([xedge[1], yedge[-1], axwidth, axheight]);
ax2.loglog(f[1:],Suk[1:-1])
ax2.loglog(f,Suk_IEC,'r')
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD (m^2/s^2/Hz)')
plt.title(fname)

# v-spectrum
ax5 =  plt.axes([xedge[1], yedge[1], axwidth, axheight]);
ax5.loglog(f[1:],Svk[1:-1])
ax5.loglog(f,Svk_IEC,'r')
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD (m^2/s^2/Hz)')

# w-spectrum
ax8 =  plt.axes([xedge[1], yedge[0], axwidth, axheight]);
ax8.loglog(f[1:],Swk[1:-1])
ax8.loglog(f,Swk_IEC,'r')
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD (m^2/s^2/Hz)')

# u-dtheta
nbins = 20
ax3 =  plt.axes([xedge[2], yedge[-1], axwidth, axheight]);
ax3.hist(dthetau,bins=nbins,normed=True)
ax3.plot(thetaPlot, fWC, 'r')
plt.xlabel('Phase Difference/pi')
plt.ylabel('uhub Phase Differences')
plt.xlim([-np.pi,np.pi])

# v-dtheta
ax6 =  plt.axes([xedge[2], yedge[1], axwidth, axheight]);
ax6.hist(dthetav,bins=nbins,normed=True)
ax6.plot(thetaPlot, fWC, 'r')
plt.xlabel('Phase Difference/pi')
plt.ylabel('vhub Phase Differences')
plt.xlim([-np.pi,np.pi])

# w-dtheta
ax9 =  plt.axes([xedge[2], yedge[0], axwidth, axheight]);
ax9.hist(dthetaw,bins=nbins,normed=True)
ax9.plot(thetaPlot, fWC, 'r')
plt.xlabel('Phase Difference/pi')
plt.ylabel('whub Phase Differences')
plt.xlim([-np.pi,np.pi])

fig.show()

if saveimg:
    plt.savefig(dname + fname + '.png')
    print fname + ' saved as image'

