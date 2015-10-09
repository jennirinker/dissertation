"""
convert airfoil geometry and aerodynamic properties in csv files into .mat 
for easier use
"""
import scipy.io as scio
import csv
import numpy as np

#AFNames = ['cylinder','s818_2702','s825_2102','s826_1602']
AFNames = ['cylinder']

for i in range(len(AFNames)):
    AFName = AFNames[i]
    
    # read in geometry
    Geo = np.empty((0,2))
    with open(AFName + '_Geo.csv','rb') as f:
        reader = csv.reader(f)
        for row in reader:
            Geo = np.vstack((Geo,
                             np.array(row,dtype=float).reshape(1,2)))
                             
    # read in lift and drag
    LiftDrag = np.empty((0,3))
    with open(AFName + '_AD.csv','rb') as f:
        reader = csv.reader(f)
        for row in reader:
            LiftDrag = np.vstack((LiftDrag,
                             np.array(row,dtype=float).reshape(1,3)))
                             
    # create dictionary and save
    AFDict = {}
    AFDict['Geometry'] = Geo
    AFDict['LiftDrag']  = LiftDrag
    scio.savemat(AFName+'.mat',AFDict)