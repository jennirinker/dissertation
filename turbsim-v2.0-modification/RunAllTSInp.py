""" Run all TurbSim input files in specified directory.

    Jenni Rinker, 08-Apr-2015
"""
import os

# directory name
dname = '2-myflag'

# get list of input files
files = os.listdir(dname)
for file in files:
    if file.endswith('.inp'):
        os.system(os.path.join(dname,'TurbSim_glin64') + ' ' + \
                  os.path.join(dname,file))
