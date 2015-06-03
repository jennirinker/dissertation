""" Print statistics to compare the metadata tables generated with matlab and 
    with python.
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import numpy as np
import scipy.io as scio
import JR_Library.data_analysis as da
import matplotlib.pyplot as plt

print_latex = 0

# =============================================================================
# load/screen the matlab metadata table

if ('cleandata_mat' not in vars()):

    # path to matlab-processed metadata table
    matpath = 'C:\\Users\\jrinker\\Dropbox\\my_publications\\' + \
        '2015-02-28_temporal coherence in data\\code\\dataStruc.mat'
    
    # load the structure
    struc = scio.loadmat(matpath)
    
    # extract the information
    dataStruc = struc['dataStruc'][0,0]
    
    # metadata table
    metadata_mat = dataStruc[0]
    
    # oddly-formatted field names
    tmp = np.squeeze(dataStruc[1])
        
    # extract readable fieldnames
    fields_mat = [str(ele[0]) for ele in tmp]
    
    # replace necessary old fieldnames with new ones for screening
    fields_mat[fields_mat.index('Cup_speed')] = 'Wind_Speed_Cup'
    fields_mat[fields_mat.index('Wind_direction')] = 'Wind_Direction'
    fields_mat[fields_mat.index('Precip_intens')] = 'Precipitation'
    fields_mat[fields_mat.index('rho_u')] = 'Concentration_u'

    # screen table
    cleandata_mat = da.screenmetadata(fields_mat,metadata_mat,'NREL')
    
    # remove nans
    nanrows = np.unique(np.nonzero(np.isnan(cleandata_mat))[0]) # rows with nans
    cleandata_mat = np.delete(cleandata_mat,nanrows,axis=0)     # delete nan rows
    
    del struc, dataStruc, tmp

# =============================================================================
# load/screen the python metadata table

if ('cleandata_py' not in vars()):

    fname = 'C:\\Users\\jrinker\\Documents\\' + \
        'GitHub\\dissertation\\data-analysis\\metadata_NREL.txt'
    fields_py, metadata_py = da.loadmetadata(fname)
    
    # screen table
    cleandata_py = da.screenmetadata(fields_py,metadata_py,'NREL')
    
    # remove nans
    nanrows = np.unique(np.nonzero(np.isnan(cleandata_py))[0]) # rows with nans
    cleandata_py = np.delete(cleandata_py,nanrows,axis=0)      # delete nan rows

# =============================================================================
# compare the height distributions

heights = [15,30,50,76,100,131]

htCol_mat = fields_mat.index('Height')
htCol_py  = fields_py.index('Height')

print 'Ht. (m)   Matlab   Python   Difference'
print '--------------------------------------'
count_mat = np.empty(len(heights))
count_py  = np.empty(len(heights))
for i in range(len(heights)):
    count_mat[i] = np.sum(cleandata_mat[:,htCol_mat] == heights[i])
    count_py[i]  = np.sum(cleandata_py[:,htCol_py] == heights[i])
    print '{0: <3}       {1:.0f}    {2:.0f}     {3:.0f}' \
        .format(heights[i],count_mat[i],count_py[i],count_mat[i]-count_py[i])
        
# =============================================================================
# compare distributions of the mean wind speed
        
UCol_mat = fields_mat.index('U')
UCol_py  = fields_py.index('Mean_Wind_Speed')        

# let's look at 15 m
height = 15
data_mat = cleandata_mat[np.where(cleandata_mat[:,htCol_mat] == height)]
data_py  = cleandata_py[np.where(cleandata_py[:,htCol_py] == height)]

# extract concentration parameters
U_mat = data_mat[:,UCol_mat]
U_py  = data_py[:,UCol_py]

# create/clear figure
plt.figure(1)
plt.clf()

# plot histograms
plt.hist(U_mat,bins=20,histtype='step',label='Matlab')
plt.hist(U_py,bins=20,histtype='step',label='Python')

# create labels, title, legend
plt.xlabel('Mean wind speed')
plt.ylabel('Histogram')
plt.title('Mean wind speed distributions for MAtlab vs. Python processing')
plt.legend()



