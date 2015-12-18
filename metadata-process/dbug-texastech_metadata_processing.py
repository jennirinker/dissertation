"""
debugging texas tech metadata processing
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import scipy.io as scio
import os
import numpy as np

# define dataset
dataset = 'texastech'

# get base directory with mat files
basedir = jr.getBasedir(dataset)

# generate and save list of mat files in directory
#listmats = jr.list_matfiles(basedir,save=1)

# pick height of interest
specs = jr.datasetSpecs(dataset)
n_t, dt, IDs  = specs['n_t'],specs['dt'],specs['IDs']      # sampling heights

# load test structure
fname = 'FT2_E05_C01_R00070_D20120121_T1010_TR.mat'
fpath = os.path.join(basedir,'2012\\01\\21',fname)
struc_hf = scio.loadmat(fpath,squeeze_me=True)

ID = 1

# manually check wind direction calculation
x  = jr.loadtimeseries(dataset,'Sonic_x',ID,struc_hf)['clean']
#y  = jr.loadtimeseries(dataset,'Sonic_y',ID,struc_hf)['clean']
#WD_sonic = np.arctan2(y,x)
#WD = sonic_offset*np.pi/180 - WD_sonic
#WDbar = np.angle(np.nanmean(np.exp(1j*WD)),deg=1)
#print('Manual wind direction [deg]: {:.1f}'.format(WDbar))
#WSbar = np.nanmean(np.sqart(x**2 + y**2))

#t = outdict['time']
#x_raw = outdict['raw']
#x_cl  = outdict['clean']
#flags = outdict['flags']

plt.figure(1)
plt.clf()

plt.plot(x,label='raw')

#plt.plot(t,x_cl,label='stored clean')
#plt.legend()

# compare calculate fieldand struc2metadata outputs
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