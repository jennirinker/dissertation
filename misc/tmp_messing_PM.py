"""
messing around with Plaine Morte data
"""
import os

baseDir = 'H:\\data\\plaine-morte_raw\\CM 2006\\Data\\03-17'
fname = 'TS_WND_2s.TOA'
fpath = os.path.join(baseDir,fname)

iLines = 100
with open(fpath,'r') as f_read:
    with open('tmp.txt','w') as f_write:
        for i in range(iLines):
            line = f_read.readline()
            f_write.write(line)