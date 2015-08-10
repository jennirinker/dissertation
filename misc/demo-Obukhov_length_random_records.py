"""
Demonstration of calculation of ``local'' Obukhov length for random records
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import scipy.io as scio
import random
import calendar, time
import numpy as np
import matplotlib.pyplot as plt

# %% =============== draw samples of records, save metadata ===================
n_recs  = 500                               # no. random records to draw
yearRng = [2012,2015]                       # range of years
fields = jr.metadataFields('NREL')          # MD fields

CSLim  = 3                          # lower cup speed limit
dir1   = 240                        # CCW edge for direction range
dir2   = 315                        # CW edge for direction range
preLim = 2.7                        # lower precipitation limit

 column indices for each value
CScol  = fields.index('Wind_Speed_Cup')
dirCol = fields.index('Wind_Direction')
preCol = fields.index('Precipitation')

i_rec = 0
metadata = np.empty((n_recs,len(fields)))
while i_rec < n_recs:
    
    # draw random time stamp/height
    year   = random.randint(yearRng[0],yearRng[1])
    month  = random.randint(1,12)
    day    = random.randint(1,31)
    hour   = random.randint(0,23)
    minute = random.randint(0,6)*10
    height = random.choice([15,30,50,76,100,131])
    time_tup = (year,month,day,hour,minute)
        
    # check if it exists
    fpath = jr.NRELtime2fpath(time_tup)
    if len(fpath) > 0:
        
        # load structure
        struc = scio.loadmat(fpath)
        
        # calculate parameters
        row = jr.extractNRELparameters(struc,height)
        
        flagCS   = row[CScol] > CSLim
        flagDir1 = row[dirCol] >= dir1
        flagDir2 = row[dirCol] <= dir2
        flagPrec = row[preCol] >= preLim
        flagNaN  = not np.isnan(row[6])
        
        if (flagCS and flagDir1 and flagDir2 and flagPrec and flagNaN):
            metadata[i_rec] = row.reshape(1,len(fields))
    
            i_rec += 1
                                    
            if not (i_rec % 5): print(i_rec)
                
# %% ============================== analysis ==================================


## get indices
rhoCol    = fields.index('Concentration_u')
zCol      = fields.index('Height')
L_intrCol = fields.index('MO_Length_interp')
L_nearCol = fields.index('MO_Length_near')
time_Col  = fields.index('Record_Time')
##
### extract values
time_UTC = metadata[:,time_Col]
time_loc = jr.timeUTC2local('NREL',time_UTC)
time_arr = jr.timeflt2arr(time_loc)
hours    = time_arr[:,3]
rho      = metadata[:,rhoCol]
z        = metadata[:,zCol]
L_intr   = metadata[:,L_intrCol]
L_near   = metadata[:,L_nearCol]
zeta     = z / L_intr

# separate by stability classes
xs_idx = np.where(zeta < -2)[0]
vs_idx = np.where(np.logical_and(zeta >= -2., zeta < -0.6))[0]
s_idx  = np.where(np.logical_and(zeta >= -0.6, zeta < -0.2))[0]
ws_idx = np.where(np.logical_and(zeta >= -0.2, zeta < -0.02))[0]
n_idx  = np.where(np.logical_and(zeta >= -0.02, zeta < 0.02))[0]
wu_idx = np.where(np.logical_and(zeta >= 0.02, zeta < 0.2))[0]
u_idx  = np.where(np.logical_and(zeta >= 0.2, zeta < 0.6))[0]
vu_idx = np.where(np.logical_and(zeta >= 0.6, zeta < 2.))[0]
xu_idx = np.where(zeta >= 2.)[0]
stab_cts = np.array([xs_idx.size,vs_idx.size,s_idx.size,ws_idx.size,n_idx.size, \
            wu_idx.size,u_idx.size,vu_idx.size,xu_idx.size])
rho_xs = rho[xs_idx]
rho_vs = rho[vs_idx]
rho_s  = rho[s_idx]
rho_ws = rho[ws_idx]
rho_n  = rho[n_idx]
rho_wu = rho[wu_idx]
rho_u  = rho[u_idx]
rho_vu = rho[vu_idx]
rho_xu = rho[xu_idx]        

# initialize figure
plt.figure(1,figsize=(11,8))
plt.clf()
plt.figtext(0.5,0.96,'500 Random Record Values from M4 Data', \
    ha='center',fontsize='x-large')
        
# PLOT 1: interpolated versus closest MO lengths
#plt.figure(1,figsize=(5,5))
ax1 = plt.axes([0.65,0.60,0.25,0.25])
plt.plot([-40000,40000],[-40000,40000],'k:')
plt.plot(L_near,L_intr,'.')
plt.xlim([-5000,5000])
plt.ylim([-5000,5000])
ax1.set_xticklabels([])
ax1.set_yticklabels([])
plt.xlabel('Closest Temp. Value')
plt.ylabel('Interp. Temp. Value')
plt.title('MO Length with Interpolated vs. Closest Temp.')

# PLOT 2: histogram of z/L
#plt.figure(2,figsize=(7,4))
#plt.clf()
ax = plt.axes([0.08,0.55,0.4,0.32])
stab_cutoffs = [-2.,-0.6,-0.2,-0.02,0.02,0.2,0.6,2]
for i in range(len(stab_cutoffs)):
    plt.plot([stab_cutoffs[i],stab_cutoffs[i]],[0,300],'k:')
plt.hist(zeta,bins=np.linspace(-4,4,31))
plt.xlabel('Local $z/L$')
plt.title('Nondimensional MO Length')
#plt.tight_layout()

# PLOT 3: histogram of stability conditions
#plt.figure(3,figsize=(7,4))
#plt.clf()
ax = plt.axes([0.08,0.10,0.4,0.32])
stab_lbls = ['XS','VS','S','WS','N','WU','U','VU','XU']
bin_cntrs = np.arange(-4,5)
cmap = plt.get_cmap('seismic')
patches = plt.bar(bin_cntrs,stab_cts,width=1.,align='center',facecolor='white')
for patch, bin_cntr in zip(patches, bin_cntrs):
    patch.set_facecolor(cmap((bin_cntr-4.)/8.+1.))
plt.xlim([-4.75,4.75])
plt.ylim([0,200])
plt.xticks(bin_cntrs, stab_lbls)
plt.title('Stability Distribution')
for count, x in zip(stab_cts, bin_cntrs):
    # Label the raw counts
    ax.annotate(str(count), xy=(x, 0), xycoords=('data', 'axes fraction'),
        xytext=(0, -18), textcoords='offset points', va='top', ha='center')

    # Label the percentages
    percent = '%0.0f%%' % (100 * float(count) / stab_cts.sum())
    ax.annotate(percent, xy=(x, 0), xycoords=('data', 'axes fraction'),
        xytext=(0, -32), textcoords='offset points', va='top', ha='center')

# PLOT 4: scaterplot of rho and L
#plt.figure(4,figsize=(7,4.5))
#plt.clf()
ax = plt.axes([0.58,0.10,0.40,0.32])
rho_mean = [np.mean(rho_xs),np.mean(rho_vs),np.mean(rho_s),np.mean(rho_ws), \
    np.mean(rho_n),np.mean(rho_wu),np.mean(rho_u),np.mean(rho_vu), \
    np.mean(rho_xu)]
rho_std  = [np.std(rho_xs),np.std(rho_vs),np.std(rho_s),np.std(rho_ws), \
    np.std(rho_n),np.std(rho_wu),np.std(rho_u),np.std(rho_vu), \
    np.std(rho_xu)]
plt.errorbar(bin_cntrs,rho_mean,yerr=rho_std,fmt='o')
plt.xlim(-4.5,4.5)
plt.xticks(bin_cntrs, stab_lbls)
plt.xlabel('Stability')
plt.ylabel('Concentration Parameter')
plt.title(r'Variation of $\rho$ with Stability')
#plt.tight_layout()
