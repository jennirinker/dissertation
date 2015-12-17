"""
debugging PM06 processing
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import scipy.io as scio
import os
import numpy as np

# define dataset
dataset = 'PM06'

# get base directory with mat files
basedir = jr.getBasedir(dataset,'G:')

# generate and save list of mat files in directory
#listmats = jr.list_matfiles(basedir,save=1)

# pick height of interest
specs = jr.datasetSpecs(dataset)
n_t, dt, IDs  = specs['n_t'],specs['dt'],specs['IDs']      # sampling heights
sonic_offset = specs['sonic_offset']

# load test structure
fname = '02_04_2006_1230_TS_WND.mat'
fpath = os.path.join(basedir,'2006\\02\\04',fname)
struc_hf = scio.loadmat(fpath)


#grp_warning = struc['grp_warning_1(1)']
#
#perc_healthy = np.sum(grp_warning != 0)/float(n_t)
#
#if (perc_healthy >= 0.95):
    
#field = 'Sonic_T'
#field = 'Humidity'
#field = 'Grp_Warning'
ID = 4

# manually check wind direction calculation
outdictx  = jr.loadtimeseries(dataset,'Sonic_x',ID,struc_hf)
x, xflags = outdictx['clean'], outdictx['flags']
outdicty  = jr.loadtimeseries(dataset,'Sonic_y',ID,struc_hf)
y, yflags = outdicty['clean'], outdicty['flags']
#WD_sonic = np.arctan2(y,x)
#WD = sonic_offset*np.pi/180 - WD_sonic
#WDbar = np.angle(np.nanmean(np.exp(1j*WD)),deg=1)
#print('Manual wind direction [deg]: {:.1f}'.format(WDbar))
WSbar = np.nanmean(np.sqrt(x**2 + y**2))

#t = outdict['time']
#x_raw = outdict['raw']
#x_cl  = outdict['clean']
#flags = outdict['flags']

#plt.figure(1)
#plt.clf()
#
#plt.plot(t,x_raw,label='raw')
#plt.plot(t,x_cl,label='stored clean')
#plt.legend()

md_fields = jr.metadataFields(dataset)
outdict = jr.calculatefield(dataset,struc_hf,ID)
row = jr.struc2metadata(dataset,struc_hf,ID)

for i_field in range(len(md_fields)):
    key = md_fields[i_field]
    print('{:20s} {:15.3f} {:15.3f}'.format(key,
                  outdict[key],
                  row[i_field]))

#print()

# load raw file for testing
#fpath_raw = 'G:\\data\\plaine-morte_raw\\CM 2006\\Data\\03-20\\TS_WND_2.TOA'
#iprint = [0,20]
#with open(fpath_raw,'r') as f:
#    for i in range(iprint[1]):
#        if (i >= iprint[0]):
#            print(f.readline())
#        else:
#            f.readline()

# calculate metadata values
#row = jr.struc2metadata(dataset,struc,height)