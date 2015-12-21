"""
Messing around with writing .bat files for Monsoon calculations
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\nwtc_python_tools'
if (libpath not in sys.path): sys.path.append(libpath)

import os
import numpy as np
import datetime
import shutil

# define turbine name, turbine dictionary
# TurbName = ['WP0.75A08V00','WP1.5A08V03','WP3.0A02V02','WP5.0A04V00']
TurbName = 'WP0.75A08V00'

SimsPerGrp = 2                              # FAST simulations per group
RunName = 'TestRun'                         # run name
DateFmt = '%a %b %d %H:%M:%S %Y'            # date formate

# directories
BatDir = os.path.join('\\\\monsoon-data\\Public\\JRinker' + \
                    '\\fast_simulations\\BatDir',RunName)
FastDir = os.path.join('\\\\monsoon-data\\Public\\JRinker\\' + \
                        'fast_simulations\\FastDir',TurbName)
TmplDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\' + \
            'nwtc_python_tools\\templates'     # loc of base templates
FastExePath = '\\\\monsoon-data\\Public\\JRinker\\' + \
                    'fast_simulations\\ExeDir\\FAST.exe'

# -----------------------------------------------------------------------------

# delete everything in BatDir if it exists
if os.path.exists(BatDir):
    shutil.rmtree(BatDir)
    
# make bat file and message directories
os.mkdir(BatDir)
MsgDir = os.path.join(BatDir,'Messages')
os.mkdir(MsgDir)

# paths to group and individual bat templates
GrpBatTmplPath = os.path.join(TmplDir,'Template_GrpBat.bat')
IndBatTmplPath = os.path.join(TmplDir,'Template_IndBat.bat')

# get list of bat files for each fast file
FastBatNames = [os.path.splitext(f)[0]for f in os.listdir(FastDir) \
            if f.endswith('.fst')]

# get number of groups
NTotSims = len(FastBatNames)
NGrps = int(np.ceil(NTotSims/float(SimsPerGrp)))
GrpDigits = int(np.ceil(np.log10(NGrps)))
IndDigits = int(np.ceil(np.log10(SimsPerGrp)))

# create RunAll file
RunAllPath = os.path.join(BatDir,'{:s}_RunAll.bat'.format(RunName))
with open(RunAllPath,'w') as RunAll:
    
    # loop through groups of FAST simulations
    iTotSims = 0
    for iGroup in range(NGrps):
        
        # define path for group bat file and message
        GrpID      = str(iGroup+1).zfill(GrpDigits)
        GrpBatName = '{:s}_Bat{:s}.bat'.format(RunName,GrpID)
        GrpBatPath = os.path.join(BatDir,GrpBatName)
        GrpBatMsgName = os.path.splitext(GrpBatName)[0] + '.mes'
        GrpBatMsgPath = os.path.join(MsgDir,GrpBatMsgName)
        
        # create time stamp line
        TimeStr   = datetime.datetime.now().strftime(DateFmt)
        TimeStamp = '::   created {:s} by J. Rinker Python code'.format(TimeStr)
        
        # open group bat files, loop through lines
        with open(GrpBatTmplPath,'r') as GrpBatTmpl:
            with open(GrpBatPath,'w') as GrpBat:
                for Line in GrpBatTmpl:
                    
                    # if line has a writeable field
                    if ('{:s' in Line):
                        
                        Field   = Line.split('{:s<')[1].split('>}')[0]
                        NewLine = ''.join(Line.split( \
                                    '<{:s}>'.format(Field))).format(eval(Field))
                        
                    else:
                        NewLine = Line
                        
                    GrpBat.write(NewLine)
                        
        
        
                iInd = 0
                # loop through simulations in that group
                while ((iTotSims < NTotSims) and (iInd < SimsPerGrp)):
                    
                    # define paths for individual bat file
                    FastName = FastBatNames[iTotSims]
                    SimID = str(iInd+1).zfill(IndDigits)
                    IndBatName = '{:s}_{:s}_S{:s}.bat'.format(RunName,FastName,SimID)
                    IndBatPath = os.path.join(BatDir,IndBatName)
                    FastOutName = os.path.splitext(IndBatName)[0] + '.out'
                    FastOutPath = os.path.join(FastDir,FastOutName)
                    
                    # ------------ write simulation to group bat file ------------
                    GrpBat.write(':: --- Simulation {:s} ---\n'.format(SimID))
                    GrpBat.write('@ECHO  Running SIM {:s} {:s}\n'.format(SimID,IndBatName))
                    GrpBat.write('@ECHO  Running SIM {:s} {:s} >> %MESSFILE%\n'.format(SimID,IndBatName))
                    GrpBat.write('@CALL "A:\\\\{:s}"\n'.format(IndBatName))
                    GrpBat.write('@ECHO    SIM {:s}  ERRORLEVEL %ERRORLEVEL% >> %MESSFILE%\n'.format(SimID))
                    GrpBat.write('\n'.format(SimID))
                    
                    # --------------- write individual bat file ---------------
                    IndBatMsgName = os.path.splitext(IndBatName)[0] + '.mes'
                    IndBatMsgPath = os.path.join(MsgDir,IndBatMsgName)
                    with open(IndBatTmplPath,'r') as IndBatTmpl:
                        with open(IndBatPath,'w') as IndBat:
                            for Line in IndBatTmpl:
                                
                                # if line has a writeable field
                                if ('{:s' in Line):
                                    
                                    Field   = Line.split('{:s<')[1].split('>}')[0]
                                    NewLine = ''.join(Line.split( \
                                                '<{:s}>'.format(Field))).format(eval(Field))
                                    
                                else:
                                    NewLine = Line
                                    
                                IndBat.write(NewLine)
                            
                    # increment simulation counters
                    iTotSims += 1
                    iInd += 1
                    
        # write group bat to run all
        Line = 'start /b job submit /scheduler:monsoon /jobname:{:s} \"\{:s}"\n'.format( \
                        GrpBatName,GrpBatPath)
        RunAll.write(Line)



