"""
Demo load time series
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import matplotlib.pyplot as plt

# Unhealthy cup wind speed data
# ***===***  6 highest values in metadata are 497,416,138,300,192,397
# ***===*** 3 of top 4 are unhealthy
time_unh = metadata[397,0]
#time_unh = (2013,01,10,23,10)
out_WS = jr.loadtimeseries('NREL','Cup_WS_10m',time_unh)
out_son = jr.loadtimeseries('NREL','Sonic_u_15m',time_unh)
t = out_WS['time']
WS_raw = out_WS['raw']
WS_cl = out_WS['clean']
u_unh = out_son['raw']

# print/plot results
print('Unhealthy Cup_WS flags: ',out_WS['flags'])

plt.figure(1)
plt.clf()

plt.subplot(311)
plt.plot(t,WS_raw)
plt.plot(t,WS_cl)
plt.xlim([-10,610])
plt.ylim([0,93])
plt.legend(['Raw Cup_WS','In-Range Cup_WS'])
plt.title('Unhealthy Wind Speed Cup Data (flags = {})'.format(out_WS['flags']))

plt.subplot(312)
plt.plot(t,u_unh)
plt.xlim([-10,610])
plt.legend(['Raw Sonic_u','Cleaned Sonic_u'])
plt.title('Unhealthy Sonic_u Data Data (flags = {})'.format(out_son['flags']))

plt.subplot(313)
plt.plot(t,WS_cl)
plt.plot(t,u_unh)
plt.title('Unhealthy Cup_WS vs. Sonic_u')
plt.legend(['Cup','Sonic_u'])

plt.tight_layout()

# healthy cup wind speed
time_hlt = metadata[1,0]
out_WS = jr.loadtimeseries('NREL','Cup_WS_10m',time_hlt)
out_son = jr.loadtimeseries('NREL','Sonic_u_15m',time_hlt)
t = out_WS['time']
WS_raw = out_WS['raw']
WS_cl = out_WS['clean']
u_unh = out_son['raw']

# print/plot results
print('Healthy Cup_WS flags: ',out_WS['flags'])

plt.figure(2)
plt.clf()

plt.subplot(311)
plt.plot(t,WS_raw)
plt.plot(t,WS_cl)
plt.xlim([-10,610])
plt.legend(['Raw Cup_WS','In-Range Cup_WS'])
plt.title('Healthy Wind Speed Cup Data (flags = {})'.format(out_WS['flags']))

plt.subplot(312)
plt.plot(t,u_unh)
plt.xlim([-10,610])
plt.legend(['Raw Sonic_u','Cleaned Sonic_u'])
plt.title('Healthy Sonic_u Data Data (flags = {})'.format(out_son['flags']))

plt.subplot(313)
plt.plot(t,WS_cl)
plt.plot(t,u_unh)
plt.title('Healthy Cup_WS vs. Sonic_u')
plt.legend(['Cup','Sonic_u'])

plt.tight_layout()