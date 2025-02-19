"""
calculating wind directions with sonic and UVW data from texas tech
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as scio

TT_sonictheta_offset = 300.3488889
TT_UVWtheta_offset = TT_sonictheta_offset - 180 - 42.96

i_t = 0
t = np.arange(30000)*0.02

h = 116
#h = 158
#h = 200

# load raw high-frequency data
fpath = 'G:\\data\\texas-tech\\2012\\01\\22\\' + \
            'FT2_E05_C01_R00098_D20120122_T0010_TR.mat'
outdict = scio.loadmat(fpath,squeeze_me=True)

# ========================= sonic data =======================================

x_s = outdict['Sonic_x_int_{:d}m'.format(h)]
y_s = outdict['Sonic_y_int_{:d}m'.format(h)]
z_z = outdict['Sonic_z_int_{:d}m'.format(h)]
x_p = outdict['UVW_x_int_{:d}m'.format(h)]
y_p = outdict['UVW_y_int_{:d}m'.format(h)]
z_p = outdict['UVW_z_int_{:d}m'.format(h)]

del outdict

WS_s = np.sqrt(x_s**2 + y_s**2)
WS_p = np.sqrt(x_p**2 + y_p**2)

plt.figure(1)
plt.clf()
plt.plot(WS_s)
plt.plot(WS_p)

#theta_loc_s = np.arctan2(y_s,x_s) * 180/np.pi
#theta_loc_p = np.arctan2(y_p,x_p) * 180/np.pi
#
#theta_c_ccw_s = 303 - 90 + theta_loc_s
#theta_c_ccw_p = 303 - 44 + theta_loc_p
#
#theta_c_s = (-theta_c_ccw_s) % 360
#theta_c_p = (-theta_c_ccw_p) % 360

#WI_s = (300 - np.arctan2(y_s,x_s) * 180/np.pi + 180.) % 360
#WI_p = (105 - np.arctan2(y_p,x_p) * 180/np.pi + 180.) % 360
#WD_s = (300 - np.arctan2(y_s,x_s) * 180/np.pi ) % 360
#WD_p = (300 + 45 - np.arctan2(y_p,x_p) * 180/np.pi ) % 360
WI_s = (300 - np.arctan2(y_s,x_s) * 180/np.pi + 180.) % 360
WI_p = (300 + 45 - np.arctan2(y_p,x_p) * 180/np.pi + 180. ) % 360


plt.figure(2)
plt.clf()
plt.plot(WI_s)
plt.plot(WI_p)
#plt.plot(WD_s)
#plt.plot(WD_p)



#u1 = dict_10mins[i_t]['Sonic_u_{:d}m'.format(h)]
#v1 = dict_10mins[i_t]['Sonic_v_{:d}m'.format(h)]
#w1 = dict_10mins[i_t]['Sonic_w_{:d}m'.format(h)]
#T1 = dict_10mins[i_t]['Sonic_T_{:d}m'.format(h)]

#vel_raw = np.concatenate((x1.reshape(x1.size,1),
#                          y1.reshape(x1.size,1),
#                          z1.reshape(x1.size,1)),axis=1)
#
## rotate time series
##vel_rot, vel_yaw, theta_rot = jr.RotateTimeSeries(vel_raw)
#
##print(np.mean(vel_raw,axis=0))
##print(np.mean(vel_yaw,axis=0))
##print(np.mean(vel_rot,axis=0))
#
## calculate different wind angles
#sonic_theta = np.arctan2(y1,x1)
#sonic_card_deg = TT_sonictheta_offset - sonic_theta*180/np.pi
#sonic_WD_deg = sonic_card_deg + 180.
#theta_avg = np.angle(np.mean(np.exp(1j*sonic_theta))) 
#
#print(theta_rot*180/np.pi,theta_avg*180/np.pi)


#plt.figure(1)
#plt.clf()

#plt.subplot(221)
#plt.plot(np.cos(sonic_theta),np.sin(sonic_theta),'.' )
#plt.xlim([-1.1,1.1])
#plt.ylim([-1.1,1.1])
#plt.title('Sonic angle in sonic axes')
#
#plt.subplot(222)
#plt.plot(np.cos(90-sonic_card_deg*np.pi/180),
#         np.sin(90-sonic_card_deg*np.pi/180),'.' )
#plt.xlim([-1.1,1.1])
#plt.ylim([-1.1,1.1])
#plt.title('Sonic angle in cardinal direction')
#plt.plot(sonic_theta*180/np.pi)

#plt.subplot(211)
#plt.plot(t,sonic_card_deg % 360)
#plt.plot(t,u1)
#plt.ylabel('Angle [deg]')
#plt.title('Sonic wind direction (cardinal)')


# ========================= UVW data =======================================

#x2 = dict_10mins[i_t]['UVW_x_int_{:d}m'.format(h)]
#y2 = dict_10mins[i_t]['UVW_y_int_{:d}m'.format(h)]
#z2 = dict_10mins[i_t]['UVW_z_int_{:d}m'.format(h)]
#u2 = dict_10mins[i_t]['UVW_u_{:d}m'.format(h)]
#v2 = dict_10mins[i_t]['UVW_v_{:d}m'.format(h)]
#w2 = dict_10mins[i_t]['UVW_w_{:d}m'.format(h)]
#T2 = dict_10mins[i_t]['Air_Temp_int_{:d}m'.format(h)]

#vel_raw = np.concatenate((x2.reshape(x2.size,1),
#                          y2.reshape(x2.size,1),
#                          z2.reshape(x2.size,1)),axis=1)
#
## convert to engineering units
##vel_raw= vel_raw * 40.26*0.447
#
## rotate time series
#vel_rot, vel_yaw, theta_rot = jr.RotateTimeSeries(vel_raw)
#
#print(np.mean(vel_raw,axis=0))
#print(np.mean(vel_yaw,axis=0))
#print(np.mean(vel_rot,axis=0))
#
## calculate different wind angles
#UVW_theta = np.arctan2(y2,x2)
#UWV_card_deg = TT_UVWtheta_offset - UVW_theta*180/np.pi
#UWV_WD_deg = UWV_card_deg + 180.
#theta_avg = np.angle(np.mean(np.exp(1j*UVW_theta))) 
#
#print(theta_rot*180/np.pi,theta_avg*180/np.pi)

#plt.subplot(212)
#plt.plot(t,sonic_card_deg % 360)
#plt.plot(t,u2)
#
#plt.subplot(212)
#plt.plot(t,T1)
#plt.plot(t,T2)

#plt.subplot(223)
#plt.plot(np.cos(UVW_theta),np.sin(UVW_theta),'.' )
#plt.xlim([-1.1,1.1])
#plt.ylim([-1.1,1.1])
#plt.title('UVW angle in UVW axes')
#
#plt.subplot(224)
##plt.ylabel('Angle [deg]')
##plt.title('Sonic wind direction (cardinal)')
#plt.plot(np.cos(90-UWV_card_deg*np.pi/180),
#         np.sin(90-UWV_card_deg*np.pi/180),'.' )
#plt.xlim([-1.1,1.1])
#plt.ylim([-1.1,1.1])
#plt.title('UVW angle in cardinal direction')

#plt.plot(t,UWV_card_deg % 360)
#plt.ylabel('Angle [deg]')
#plt.title('UVW wind direction (cardinal)')
#plt.plot(UVW_theta*180/np.pi)

#plt.subplot(325)
#plt.plot(t, sonic_card_deg - UWV_card_deg)
#plt.title('Difference in sonic and UVW angles')
#plt.xlabel('Time [s]')
#plt.ylabel('Angle [deg]')

#print(np.mean(sonic_card_deg - UWV_card_deg))

#print(np.mean(T2/T1))
#print(np.mean(T2-T1))
#
#plt.tight_layout()
