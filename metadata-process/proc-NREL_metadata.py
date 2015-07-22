"""
Script to do the following:
    1) iterate through NREL 20-Hz .mat files;
    2) load and clean time series;
    3) save screened/cleaned data in metadata table;
    4) save processed metadata.
"""
import scipy.io as scio
import time

start = time.time()

fpath = 'G:\\data\\nrel-20Hz\\2012\\02\\13\\02_13_2012_17_30_00_038.mat'
struc = scio.loadmat(fpath)
height = 15

heights = np.array([15,30,50,76,100,131])   # sampling heights
CSLim  = 3                          # lower cup speed limit
dir1   = 240                        # CCW edge for direction range
dir2   = 315                        # CW edge for direction range
preLim = 2.7                        # lower precipitation limit
fields = jr.metadataFields('NREL')     # list of metadata columns

# column indices for each value
CScol  = fields.index('Wind_Speed_Cup')
dirCol = fields.index('Wind_Direction')
preCol = fields.index('Precipitation')

row = jr.extractNRELparameters(struc,height)


flagCS   = row[CScol] > CSLim
flagDir1 = row[dirCol] >= dir1
flagDir2 = row[dirCol] <= dir2
flagPrec = row[preCol] >= preLim
flagNaN  = not np.isnan(row[6])

print(flagCS , flagDir1 , flagDir2 , flagPrec , flagNaN)

tmp = np.empty([1,len(fields)])
tmp[0,:] = row.reshape(1,len(fields))