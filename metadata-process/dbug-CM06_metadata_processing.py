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
dataset = 'CM06'

# get base directory with mat files
basedir = jr.getBasedir(dataset)

# generate and save list of mat files in directory
#listmats = jr.list_matfiles(basedir,save=1)

# pick height of interest
n_t, dt, IDs  = jr.datasetSpecs(dataset)      # sampling heights

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

#outdict  = jr.loadtimeseries(dataset,field,ID,struc_hf)
#t = outdict['time']
#x_raw = outdict['raw']
#x_cl  = outdict['clean']
#flags = outdict['flags']
#
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