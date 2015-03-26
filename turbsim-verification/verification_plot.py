""" Script to load a turbsim file and create a plot
	 with subplots for verification. """

import pyts.io.main as io
import matplotlib.pyplot as plt
import numpy as np
import jr.TS_Verification_Library as jr

# name of file to load
dname   = '2-periodic/TS/'
fname   = 'IEC_scale_10pts';

# save image in directory?
saveimg = 1

# construct total file path
inpname = dname + fname + '.inp';
outname = dname + fname + '.wnd';

# read file
tsout = io.readModel(outname);
tsin  = jr.readInput(inpname);

# useful values
y = tsout.grid.y;               # y-grid vector
# z = tsout.grid.z[::-1];         # reverse to match TS output
z = tsout.grid.z;         # reverse to match TS output  
[Y, Z] = np.meshgrid(y, z);     # grid arrays
n_t = tsout.uhub.size;          # (grid.n_t is not correct)
n_f = np.ceil((n_t-1)/2+1);     # unique components counting DC
rsep = min(tsout.grid.dy,\
           tsout.grid.dz);      # coherence sep distance

# ============== Data ==============

# calculate velocity and turb intesn profiles
U  = jr.velocityProfile(tsout);
Ti = tsout.Ti;
		
# calculate spatial coherence
[f,Coh] = jr.calculateTurbSimSC(outname,rsep);

# calculate hub-height spectra
[Su,Sv,Sw] = jr.hubHeightPSDs(tsout);

# phase differences
[dthetau,dthetav,dthetaw] = \
    jr.hubHeightPhaseDiffs(tsout);

# ============== Theory ==============
U_IEC  = jr.IEC_VelProfile(z,tsout);
TI_IEC = jr.IEC_TiProfile(z,tsin,tsout);
Coh_IEC = jr.IEC_SpatialCoherence(rsep,f,tsout);
[Su_IEC,Sv_IEC,Sw_IEC] = jr.IEC_PSDs(f,tsin,tsout);

# ============== Plot ==============
axwidth = 0.23;
axheight = 0.23;
xedge = [0.08, 0.41, 0.74];
yedge = [0.09, 0.41, 0.73];

fig = plt.figure(1,figsize=(9.05,7.14));
fig.clf();
plt.ion()

# velocity profile
ax1 = plt.axes([xedge[0], yedge[-1], axwidth, axheight])
ax1.scatter( U, Z )
ax1.plot( U_IEC, z, 'r')
plt.xlabel('Velocity (m/s)')
plt.ylabel('Height (m)')

# turbulence intensity profile
ax4 = plt.axes([xedge[0], yedge[1], axwidth, axheight]);
ax4.scatter( Ti, Z )
ax4.plot( TI_IEC, z, 'r')
plt.xlabel('Turb. Intens (-)')
plt.ylabel('Height (m)')

# spatial coherence
ax7 = plt.axes([xedge[0], yedge[0], axwidth, axheight]);
ax7.semilogx(f,Coh);
ax7.semilogx(f,Coh_IEC,'r')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Spatial Coherence (-)')

# u-spectrum
ax2 = plt.axes([xedge[1], yedge[-1], axwidth, axheight]);
ax2.loglog(f,Su)
ax2.loglog(f,Su_IEC,'r')
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD (m^2/s^2/Hz)')

# v-spectrum
ax5 =  plt.axes([xedge[1], yedge[1], axwidth, axheight]);
ax5.loglog(f,Sv)
ax5.loglog(f,Sv_IEC,'r')
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD (m^2/s^2/Hz)')

# w-spectrum
ax8 =  plt.axes([xedge[1], yedge[0], axwidth, axheight]);
ax8.loglog(f,Sw)
ax8.loglog(f,Sw_IEC,'r')
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD (m^2/s^2/Hz)')

# u-dtheta
ax3 =  plt.axes([xedge[2], yedge[-1], axwidth, axheight]);
ax3.hist( dthetau / np.pi, 20 );
plt.xlabel('Phase Difference/pi')
plt.ylabel('uhub Phase Differences')

# v-dtheta
ax6 =  plt.axes([xedge[2], yedge[1], axwidth, axheight]);
ax6.hist( dthetav / np.pi, 20 );
plt.xlabel('Phase Difference/pi')
plt.ylabel('vhub Phase Differences')

# w-dtheta
ax9 =  plt.axes([xedge[2], yedge[0], axwidth, axheight]);
ax9.hist( dthetaw / np.pi, 20 );
plt.xlabel('Phase Difference/pi')
plt.ylabel('whub Phase Differences')

fig.show()

if saveimg:
    plt.savefig(dname + fname + '.png')

