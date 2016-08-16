"""
Scatterplot sigu vs U for 1 m and 200 m
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)
import JR_Library.main as jr

import numpy as np
import os, json
import matplotlib.pyplot as plt


# choose which dataset
dataset = 'texastech'

# define directory where wind parameters are stored (unused for matlab)
basedir = 'C:\\Users\\jrinker\\Dropbox\\research\\' + \
                  'processed_data'

FigNum = 1
FigSize = (6.0,3.0)

# plot style
plt.style.use(jr.stylepath('duke_paper'))

colors = ['#235F9C', '#C0504D', '#F79646', '#8064A2', \
            '#9BBB59', '#8C564B', '#17BECF', '#BCBD22', '#7F7F7F', '#E377C2', \
            '#262626', '#FFD960']

# -----------------------------------------------------------------------------

# ---------- load data ----------
dist_type = 'comp'

if ('mat' in dataset):
    fields, clean = jr.loadNRELmatlab()
    dataset_flag = 'NREL'
else:
    fields, clean = jr.loadmetadata(dataset)
    dataset_flag = dataset

# screen metadata, get measurement heights and columns for data
screen = jr.screenmetadata(fields,clean,dataset_flag)
IDs     = jr.datasetSpecs(dataset_flag)['IDs']
htCol   = fields.index('ID')
UCol    = fields.index('Mean_Wind_Speed')
siguCol = fields.index('Sigma_u')

ID_idcs = [0,-1]

fig = plt.figure(1,figsize=FigSize)
plt.clf()

for iHt in range(len(ID_idcs)):
    ID_idx = ID_idcs[iHt]
    
    idcs = screen[:,htCol] == IDs[ID_idx]
    U = screen[idcs,UCol]
    sigu = screen[idcs,siguCol]
    
    plt.subplot(len(ID_idcs),1,iHt+1)
    plt.scatter(U,sigu)
