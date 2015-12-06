"""
calculating wind directions with sonic and UVW data from texas tech
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
import matplotlib.pyplot as plt

TT_sonictheta_offset = 300.3488889
TT_UVWtheta_offset = TT_sonictheta_offset - 180 - 42.96

i_t = 0
t = np.arange(30000)*0.02

# ========================= sonic data =======================================

x1 = dict_10mins[i_t]['Sonic_x_int_200m']
y1 = dict_10mins[i_t]['Sonic_y_int_200m']
z1 = dict_10mins[i_t]['Sonic_z_int_200m']

vel_raw = np.concatenate((x1.reshape(x1.size,1),
                          y1.reshape(x1.size,1),
                          z1.reshape(x1.size,1)),axis=1)

# rotate time series
vel_rot, vel_yaw = jr.RotateTimeSeries(vel_raw)

#print(np.mean(vel_raw,axis=0))
#print(np.mean(vel_yaw,axis=0))
#print(np.mean(vel_rot,axis=0))

# calculate different wind angles
sonic_theta = np.arctan2(y1,x1)
sonic_card_deg = TT_sonictheta_offset - sonic_theta*180/np.pi
sonic_WD_deg = sonic_card_deg + 180.
theta_rot = np.arctan2(np.mean(y1),np.mean(x1))
theta_avg = np.angle(np.mean(np.exp(1j*sonic_theta))) 

print(theta_rot*180/np.pi,theta_avg*180/np.pi)


print(WD_avg )

plt.figure(1)
plt.clf()

plt.subplot(221)
plt.plot(np.cos(sonic_theta),np.sin(sonic_theta),'.' )
plt.xlim([-1.1,1.1])
plt.ylim([-1.1,1.1])
plt.title('Sonic angle in sonic axes')

plt.subplot(222)
#plt.plot(t,sonic_card_deg % 360)
#plt.ylabel('Angle [deg]')
#plt.title('Sonic wind direction (cardinal)')
plt.plot(np.cos(90-sonic_card_deg*np.pi/180),
         np.sin(90-sonic_card_deg*np.pi/180),'.' )
plt.xlim([-1.1,1.1])
plt.ylim([-1.1,1.1])
plt.title('Sonic angle in cardinal direction')
#plt.plot(sonic_theta*180/np.pi)
#plt.plot(vel_rot[:,0])

# ========================= UVW data =======================================

x2 = dict_10mins[i_t]['UVW_x_int_200m']
y2 = dict_10mins[i_t]['UVW_y_int_200m']
z2 = dict_10mins[i_t]['UVW_z_int_200m']

vel_raw = np.concatenate((x2.reshape(x2.size,1),
                          y2.reshape(x2.size,1),
                          z2.reshape(x2.size,1)),axis=1)

# convert to engineering units
vel_raw= vel_raw * 40.26*0.447

# rotate time series
vel_rot, vel_yaw = jr.RotateTimeSeries(vel_raw)

print(np.mean(vel_raw,axis=0))
print(np.mean(vel_yaw,axis=0))
print(np.mean(vel_rot,axis=0))

# calculate different wind angles
UVW_theta = np.arctan2(y2,x2)
UWV_card_deg = TT_UVWtheta_offset - UVW_theta*180/np.pi
UWV_WD_deg = UWV_WD_deg + 180.
theta_rot = np.arctan2(np.mean(y2),np.mean(x2))
theta_avg = np.angle(np.mean(np.exp(1j*UVW_theta))) 

print(theta_rot*180/np.pi,theta_avg*180/np.pi)


print(WD_avg )

plt.subplot(223)
plt.plot(np.cos(UVW_theta),np.sin(UVW_theta),'.' )
plt.xlim([-1.1,1.1])
plt.ylim([-1.1,1.1])
plt.title('UVW angle in UVW axes')

plt.subplot(224)
#plt.plot(t,sonic_card_deg % 360)
#plt.ylabel('Angle [deg]')
#plt.title('Sonic wind direction (cardinal)')
plt.plot(np.cos(90-UWV_card_deg*np.pi/180),
         np.sin(90-UWV_card_deg*np.pi/180),'.' )
plt.xlim([-1.1,1.1])
plt.ylim([-1.1,1.1])
plt.title('UVW angle in cardinal direction')

#plt.plot(t,UWV_card_deg % 360)
#plt.ylabel('Angle [deg]')
#plt.title('UVW wind direction (cardinal)')
#plt.plot(UVW_theta*180/np.pi)
#plt.plot(vel_rot[:,0])

#plt.subplot(325)
#plt.plot(t, sonic_card_deg - UWV_card_deg)
#plt.title('Difference in sonic and UVW angles')
#plt.xlabel('Time [s]')
#plt.ylabel('Angle [deg]')

print(np.mean(sonic_card_deg - UWV_card_deg))

plt.tight_layout()
