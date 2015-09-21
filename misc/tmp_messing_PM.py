"""
messing around with Plaine Morte data
"""
import os

baseDir = 'H:\\data\\plaine-morte_raw\\CM 2008\\Data\\01 04 2008'
fname = 'TS_WND_1.TOA'
fpath = os.path.join(baseDir,fname)

iLines = 10
with open(fpath,'r') as f_read:
    with open('tmp.txt','w') as f_write:
        for i in range(iLines):
            line = f_read.readline()
            f_write.write(line)