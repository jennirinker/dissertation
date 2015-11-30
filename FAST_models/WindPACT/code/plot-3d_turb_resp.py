"""
3d turbine visualization
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import matplotlib.pyplot as plt
import os
import pyts.io.main as io
import json
import numpy as np

# set directory and turbine name
#turb_dir,turb_name = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP0.75A08V00','WP0.75A08V00'
#turb_dir,turb_name = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP0.75A08V00_newGBR','WP0.75A08V00'
#turb_dir,turb_name = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP1500_FAST_v7','WP1500'
turb_dir,turb_name = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
    'dissertation\\FAST_models\\FAST7\\WP1.5A08V03','WP1.5A08V03'
#turb_dir,turb_name = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP3.0A02V02','WP3.0A02V02'
#turb_dir,turb_name = 'C:\\Users\\jrinker\\Documents\\GitHub\\' + \
#    'dissertation\\FAST_models\\FAST7\\WP5.0A04V00','WP5.0A04V00'
#fileID = '00000'
#turb_dir,turb_name = 'C:\\Users\\jrinker\\Documents\\GitHub' + \
#        '\\dissertation\\FAST_models\\verification\\WP0.75A08V00','WP0.75A08V00'
#fileID = ['24134', '24142', '42331'][0]
fileID = ['24211', '24214', '24220'][0]

wind_dir = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation\\' + \
            'FAST_models\\wind_files\\WP0.75A08V00'
wind_fname = 'WP0.75A08V00_24134.bts'

fignum = 2

t_lim = [240,300]

# load turbine dictionary
fTDicturb_name = os.path.join(turb_dir,'parameters',turb_name+'_Dict.dat')
with open(fTDicturb_name,'r') as f:
    turb_dict = json.load(f)

# load input file
wind_fpath = os.path.join(wind_dir,wind_fname)
tsout = io.readModel(wind_fpath)            # read file
u = tsout.u[::-1,:,:]
t_ts = tsout.time

# load fast file
FASTfname = turb_name+'_'+fileID
FASTfpath = os.path.join(turb_dir,FASTfname+'.out')
FAST = jr.ReadFASTFile(FASTfpath)
Fields,Data = FAST['Fields'],FAST['Data']
t = Data[:,FAST['Fields'].index('Time')]

# initialize figure
fig = plt.figure(fignum,figsize=(10.,6.))
#fig.clf()


# intermediate paramters (time independent)
i_fs_start = np.abs(t-t_lim[0]).argmin()
i_fs_end   = np.abs(t-t_lim[1]).argmin()
i_ts_start = np.abs(t_ts-t_lim[0]).argmin()
i_ts_end   = np.abs(t_ts-t_lim[1]).argmin()

u_plots = u[:,:,i_ts_start:i_ts_end]
levels  = np.linspace(np.floor(u_plots.min()),np.ceil(u_plots.max()),31)

y,z = tsout.grid.y, tsout.grid.z
Y,Z = np.meshgrid(y,z)
y_ticklist = jr.grid2ticklist(y)
z_ticklist = jr.grid2ticklist(z)

t_plot = t[i_fs_start:i_fs_end]
Azimuth = Data[i_fs_start:i_fs_end,Fields.index('Azimuth')] + 90
OoPDefl = np.empty((Azimuth.size,3))
IPDefl  = np.empty((Azimuth.size,3))
OoPDefl[:,0] = Data[i_fs_start:i_fs_end,Fields.index('OoPDefl1')]
OoPDefl[:,1] = Data[i_fs_start:i_fs_end,Fields.index('OoPDefl2')]
OoPDefl[:,2] = Data[i_fs_start:i_fs_end,Fields.index('OoPDefl3')]
IPDefl[:,0] = Data[i_fs_start:i_fs_end,Fields.index('IPDefl1')]
IPDefl[:,1] = Data[i_fs_start:i_fs_end,Fields.index('IPDefl2')]
IPDefl[:,2] = Data[i_fs_start:i_fs_end,Fields.index('IPDefl3')]
TTDspFA = Data[i_fs_start:i_fs_end,Fields.index('TTDspFA')]
TTDspSS = Data[i_fs_start:i_fs_end,Fields.index('TTDspSS')]

max_defl = np.ceil(OoPDefl.max())
min_defl = np.floor(OoPDefl.min())

HubHt = turb_dict['Nacelle']['HubHeight']
HubRd = turb_dict['Rotor']['HubDiam']/2.
RotRd = turb_dict['Rotor']['RotDiam']/2.
BldLn = RotRd - HubRd

turb_ylim = [0,np.ceil((HubHt+RotRd)/10.)*10]

angles = np.linspace(0,np.pi*2,300)
bld_offsts = [0,2*np.pi/3,4*np.pi/3]

IP_scale = 1
OP_scale = 1
colors = ['b','g','r']

# loop through FAST time steps 
#for i_t in range(Azimuth.size):
i_t = 0
    
print('Plotting {:d} of {:d}'.format(i_t,Azimuth.size))

fig.clf()
ax_t = plt.axes([0.07,0.5,0.13,0.2])    # turbulence axes
ax_f = plt.axes([0.3,0.27,0.3,0.7])     # front-facing turbine
ax_l = plt.axes([0.65,0.27,0.3,0.7])    # left-facing turbine
ax_b = plt.axes([0.07,0.05,0.88,0.15])  # blade response

RotAng = Azimuth[i_t]/180.*np.pi

# plot stuff
fig.text(0.02,0.95,'t = {:.1f} s'.format(t_plot[i_t]),fontsize='large')

# =============== turbulent field ===============
i_ts = np.abs(t_ts-t_plot[i_t]).argmin()
cnt = ax_t.contourf(Y,Z,u[:,:,i_ts],
              cmap='RdBu_r',levels=levels)
plt.colorbar(cnt,ax=ax_t)
ax_t.set_title('Wind input',fontsize='medium')
ax_t.set_xticks([tup[0] for tup in y_ticklist])
ax_t.set_xticklabels([tup[1] for tup in y_ticklist])
ax_t.set_yticks([tup[0] for tup in z_ticklist])
ax_t.set_yticklabels([tup[1] for tup in z_ticklist])

# =============== blade deflections ===============
ax_b.plot([t_plot[i_t],t_plot[i_t]],[min_defl,max_defl],'k:')
ax_b.plot(t_plot,OoPDefl)
ax_b.set_ylim([min_defl,max_defl])
ax_b.set_ylabel('Blade OoP Defl')

# =============== front-facing turbine ===============

# undeflected tower and rotor
ax_f.plot([0,0],[0,HubHt],'k:')
ax_f.plot(RotRd*np.cos(angles),RotRd*np.sin(angles)+HubHt,'k:')
ax_f.plot(HubRd*np.cos(angles),
          HubRd*np.sin(angles)+HubHt,
          'k:')

# deflected tower and hub
ax_f.plot(HubRd*np.cos(angles)-TTDspSS[i_t],
          HubRd*np.sin(angles)+HubHt,
          'k')
ax_f.plot([0,-TTDspSS[i_t]],[0,HubHt-HubRd],'k')

# undeflected blades
for i_bl in range(3):
    theta_und = RotAng + bld_offsts[i_bl]
    x0 = TTDspSS[i_t] + HubRd*np.cos(theta_und)
    y0 = HubHt        + HubRd*np.sin(theta_und)
    xf = x0 + BldLn*np.cos(theta_und)
    yf = y0 + BldLn*np.sin(theta_und)
    ax_f.plot([x0,xf],[y0,yf],colors[i_bl]+'--')

# deflected blades
for i_bl in range(3):
    theta_und = RotAng + bld_offsts[i_bl]
    x0 = TTDspSS[i_t] + HubRd*np.cos(theta_und)
    y0 = HubHt        + HubRd*np.sin(theta_und)
    xf_und = x0 + BldLn*np.cos(theta_und)
    yf_und = y0 + BldLn*np.sin(theta_und)
    xf = xf_und + IP_scale*IPDefl[i_t,i_bl]*np.sin(theta_und)
    yf = yf_und - IP_scale*IPDefl[i_t,i_bl]*np.cos(theta_und)
    ax_f.plot([x0,xf],[y0,yf],colors[i_bl])

ax_f.set_ylim(turb_ylim)
ax_f.text(0.105, 0.1,r'$\vec{u}$',
          ha='center',transform=ax_f.transAxes)
ax_f.text(0.1, 0.06,r'$\times$',
          ha='center',fontsize='large', transform=ax_f.transAxes)

# =============== left-facing turbine ===============

# make positive deflection to the left

# undeflected tower and rotor
ax_l.plot([0,0],[0,HubHt],'k:')
ax_l.plot([min_defl,max_defl],[HubHt-RotRd,HubHt-RotRd],'k:')
ax_l.plot([min_defl,max_defl],[HubHt+RotRd,HubHt+RotRd],'k:')

# deflected tower and hub
ax_l.plot([0,TTDspFA[i_t]],[0,HubHt],'k')

# undeflected blades
for i_bl in range(3):
    theta_und = RotAng + bld_offsts[i_bl]
    x0 = TTDspFA[i_t]
    y0 = HubHt        + HubRd*np.sin(theta_und)
    xf = x0 
    yf = y0 + BldLn*np.sin(theta_und)
    ax_l.plot([x0,xf],[y0,yf],colors[i_bl]+'--')

# deflected blades
for i_bl in range(3):
    theta_und = RotAng + bld_offsts[i_bl]
    x0 = TTDspFA[i_t]
    y0 = HubHt        + HubRd*np.sin(theta_und)
    xf_und = x0 
    yf_und = y0 + BldLn*np.sin(theta_und)
    xf     = xf_und + OP_scale*OoPDefl[i_t,i_bl]
    yf     = yf_und
    ax_l.plot([x0,xf],[y0,yf],colors[i_bl])

ax_l.set_ylim(turb_ylim)
ax_l.set_xlim([min_defl,max_defl])
ax_l.text(0.105, 0.1,r'$\vec{u}$',
          ha='center',transform=ax_l.transAxes)
ax_l.annotate('', xy=(0.18, 0.07), xycoords='axes fraction',
                xytext=(0.05, 0.07), textcoords='axes fraction',
                arrowprops=dict(arrowstyle="->")
                )

# save plot to stich movie
#plt_name = '{:s}_{:d}'.format(turb_name,i_t)
#fig.savefig(os.path.join(turb_dir,'movie',plt_name)+'.png')


