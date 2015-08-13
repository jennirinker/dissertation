"""
testing fluela data
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
    
import JR_Library.main as jr
import numpy as np

dataset = 'fluela'

basedir = jr.getBasedir(dataset)

timestamp = (2009,12,21,10,0)

fpath = jr.time2fpath(dataset,timestamp)

struc = scio.loadmat(fpath)

heights = [36,54,75]
dt, n_t = 0.1, 6000
t = np.arange(n_t)*dt
df = 1./(dt*n_t)
f = np.arange(jr.uniqueComponents(n_t))*df

datfield = jr.field2datfield(dataset,'Sonic_Cup',36)
u_cup = struc[datfield][0,0][0]
print(np.nanmean(u_cup))

plt.figure(1)
plt.clf()
field = 'Sonic_T'
for i in range(len(heights)):
    ht = heights[i]
    outdict = jr.loadtimeseries(dataset,field,ht,timestamp)
    t = outdict['time']
    u_cln = outdict['clean']
    tau = jr.calculateKaimal(u_cln,dt)
    U = np.fft.rfft(u_cln)/n_t
    sig = np.nanstd(u_cln)
    S_kai = jr.KaimalSpectrum(f,tau,sig)
    S_dat = np.power(np.abs(U),2)*2/df
    
    plt.subplot(3,1,i+1)
#    plt.loglog(f,S_dat,'.')
#    plt.loglog(f,S_kai)
    plt.plot(t,u_cln)
