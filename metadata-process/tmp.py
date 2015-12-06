"""
calculating wind directions with sonic and UVW data from texas tech
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import numpy as np
import matplotlib.pyplot as plt

TT_WD_offset = -60

i_t = 0
    

# ========================= sonic data =======================================

x = dict_10mins[i_t]['Sonic_x_int_4m']
y = dict_10mins[i_t]['Sonic_y_int_4m']
z = dict_10mins[i_t]['Sonic_z_int_4m']

vel_raw = np.concatenate((x.reshape(x.size,1),
                          y.reshape(x.size,1),
                          z.reshape(x.size,1)),axis=1)

# rotate time series
vel_rot, vel_yaw = jr.RotateTimeSeries(vel_raw)

print(np.mean(vel_raw,axis=0))
print(np.mean(vel_yaw,axis=0))
print(np.mean(vel_rot,axis=0))

sonic_theta = np.arctan2(y,x)

# calculate different wind angles
theta_rot = np.arctan2(np.mean(y),np.mean(x))
theta_avg = np.angle(np.mean(np.exp(1j*sonic_theta))) 

print(theta_rot*180/np.pi,theta_avg*180/np.pi)

# calculate average wind direction (in terms of cardinal directions)
WD_avg = 360 + TT_WD_offset - theta_avg*180/np.pi

print(WD_avg)

plt.figure(1)
plt.clf()
plt.plot(vel_rot[:,0])

# ========================= UVW data =======================================

x = dict_10mins[i_t]['UVW_x_int_4m']
y = dict_10mins[i_t]['UVW_y_int_4m']
z = dict_10mins[i_t]['UVW_z_int_4m']

vel_raw = np.concatenate((x.reshape(x.size,1),
                          y.reshape(x.size,1),
                          z.reshape(x.size,1)),axis=1)

# rotate time series
vel_rot, vel_yaw = jr.RotateTimeSeries(vel_raw)

print(np.mean(vel_raw,axis=0))
print(np.mean(vel_yaw,axis=0))
print(np.mean(vel_rot,axis=0))




