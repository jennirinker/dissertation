"""
Write bat files for monsoon simulation: just TSRun and FSRun -- split into groups
"""
import sys
#libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\nwtc_python_tools'
#if (libpath not in sys.path): sys.path.append(libpath)
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import os
import datetime
import shutil
import numpy as np

# define turbine name, turbine dictionary
TurbNames = ['WP0.75A08V00','WP1.5A08V03','WP3.0A02V02','WP5.0A04V00']
#TurbNames = ['WP0.75A08V00']


# run specifics
DateFmt = '%a %b %d %H:%M:%S %Y'            # date format
#RunName = 'TestBig'                        # run name
RunName = 'BigRun1'                        # run name
RunParms = jr.RunName2WindParms(RunName)
URefs,Is,Ls,rhos,n_dups = RunParms['URefs'],RunParms['Is'],RunParms['Ls'], \
                            RunParms['rhos'],RunParms['n_dups']

# simulations per group
SimsPerGrp = 20

# directories
BaseBatDir = os.path.join('\\\\monsoon-data\\Public\\JRinker' + \
                    '\\fast_simulations\\BatDir',RunName)
TmplDir = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\' + \
            'nwtc_python_tools\\templates'     # loc of base templates
ExeDir = '\\\\monsoon-data\\Public\\JRinker\\' + \
                    'fast_simulations\\ExeDir'
BaseWindDir = os.path.join('\\\\monsoon-data\\Public\\JRinker' + \
                    '\\fast_simulations\\WindDir',RunName)
BaseFastDir = os.path.join('\\\\monsoon-data\\Public\\JRinker\\' + \
            'fast_simulations\\FastDir',RunName)
BaseModlDir = '\\\\monsoon-data\\Public\\JRinker\\fast_simulations\\ModlDir'
ExeNames = ['TurbSim.exe','FAST.exe']
ExeIDs   = ['RunTS','RunFS']
                    
# -----------------------------------------------------------------------------

print('\nDeleting base directories...\n')

# clear BaseBatDir
if os.path.exists(BaseBatDir): shutil.rmtree(BaseBatDir)
os.mkdir(BaseBatDir)
TempDir = os.path.join(BaseBatDir,'tmp')
os.mkdir(TempDir)
if os.path.exists(BaseFastDir): shutil.rmtree(BaseFastDir)
os.mkdir(BaseFastDir)
    
# paths to bat template
BatTmplPath = os.path.join(TmplDir,'Template.bat')

# number of groups
NumSims = len(URefs)*len(Is)*len(Ls)*len(rhos)*n_dups
NumGrps = int(np.ceil(NumSims/float(SimsPerGrp)))

# get file IDs
FileIDs = []
for iU,iI,iL,iR in [(a,b,c,d) for a in range(len(URefs)) \
                            for b in range(len(Is)) \
                            for c in range(len(Ls)) \
                            for d in range(len(rhos))]:
                                
    # loop through duplicates               
    for iS in range(n_dups):

        # create file ID
        FileID = format(iU,'x')+format(iI,'x')+ \
                    format(iL,'x')+format(iR,'x')+ \
                    format(iS,'x') 
                    
        FileIDs.append(FileID)

# loop through RunTS and RunFS
for iExe in range(len(ExeNames)):
    
    ExeName = ExeNames[iExe]
    ExeID   = ExeIDs[iExe]
    ExePath = os.path.join(ExeDir,ExeName)
    TimeStr   = datetime.datetime.now().strftime(DateFmt)
    BatCmnt = '::   bat file to run {:s} '.format(ExeName) + \
                 'created {:s} by J. Rinker Python code'.format(TimeStr)

    print('Writing bat files for \"{:s}\"'.format(ExeID))
    
    # loop through turbines
    GrpID = 0
    for TurbName in TurbNames:
        
        # define folder for wind and FAST files
        BatDir  = os.path.join(BaseBatDir,TurbName)
        WindDir = os.path.join(BaseWindDir,TurbName)
        FastDir = os.path.join(BaseFastDir,TurbName)
        ModlDir = os.path.join(BaseModlDir,TurbName)
    
        # create turbine-specific directories if they don't exist
        if not os.path.exists(BatDir): os.mkdir(BatDir)
        if not os.path.exists(WindDir): os.mkdir(WindDir)
        if not os.path.exists(FastDir): os.mkdir(FastDir)
        MsgDir = os.path.join(BatDir,'Messages')
        if not os.path.exists(MsgDir): os.mkdir(MsgDir)
            
        # copy pitch file to FastDir
        PitchStart = os.path.join(ModlDir,'pitch.ipt')
        PitchEnd   = os.path.join(FastDir,'pitch.ipt')
        shutil.copyfile(PitchStart,PitchEnd)
                    
                    
        # loop through group bat files
        for iGrp in range(NumGrps):
            
            # get starting and ending indices from all simulations
            iStart = iGrp * SimsPerGrp
            iEnd   = min((iGrp+1)*SimsPerGrp,NumSims)
            
            # open RunAll file
            RunAllName = '{:s}_{:s}_Group{:d}.bat'.format(RunName,ExeID,GrpID)
            RunAllPath = os.path.join(BaseBatDir,RunAllName)
                                  
            with open(RunAllPath,'w') as RunAll:
                
                # unmap used drives if first group
                if iGrp == 0:
                    RunAll.write('net use A: /delete\n')
                    RunAll.write('net use S: /delete\n')
                RunAll.write('@ECHO Off\n')
                
                # loop through simulations
                for iSim in range(iStart,iEnd):
                    
                    # get file ID and turbine name
                    FileID   = FileIDs[iSim]
                    
                    # set directories
                    BatDir  = os.path.join(BaseBatDir,TurbName)
                    WindDir = os.path.join(BaseWindDir,TurbName)
                    FastDir = os.path.join(BaseFastDir,TurbName)
                    
                    # bat and message file names and paths
                    BatName    = '{:s}_{:s}_{:s}.bat'.format(TurbName,FileID,ExeID)
                    BatMsgName = '{:s}_{:s}_{:s}.mes'.format(TurbName,FileID,ExeID)
                    BatPath    = os.path.join(BatDir,BatName)
                    
                    # paths to TurbSim and FAST input and output files
                    TurbSimInpName = '{:s}_{:s}.inp'.format(TurbName,FileID)
                    TurbSimOutName = '{:s}_{:s}.bts'.format(TurbName,FileID)
                    FastInpName    = '{:s}_{:s}.fst'.format(TurbName,FileID)
                    TempFileName   = 'Done{:s}.tmp'.format(FileID)
                    TurbSimInpPath = os.path.join(WindDir,TurbSimInpName)
                    TurbSimOutPath = os.path.join(WindDir,TurbSimOutName)
                    FastInpPath    = os.path.join(FastDir,FastInpName)
                    TempFilePath   = os.path.join(TempDir,TempFileName)
        
                    # exectuable command
                    if (iExe == 0):
                        ExeCmnd = '%EXE% {:s}'.format(TurbSimInpPath)
                    elif (iExe == 1):
                        ExeCmnd = '%EXE% {:s}'.format(FastInpPath)
        
                    # open the full run bat file, loop through lines
                    with open(BatTmplPath,'r') as Tmpl:
                        with open(BatPath,'w') as Bat:
                            for Line in Tmpl:
                                
                                # if line has a writeable field
                                if ('{:s' in Line):
                                    Field   = Line.split('{:s<')[1].split('>}')[0]
                                    NewLine = ''.join(Line.split( \
                                                '<{:s}>'.format(Field))).format(eval(Field))
                                    
                                else:
                                    NewLine = Line
                                    
                                Bat.write(NewLine)
                                
                                
                            
                    # write group bat to run all
                    Line = 'start /b job submit /scheduler:monsoon' + \
                        ' /jobname:{:s} \"{:s}"\n'.format(BatName,BatPath)
                    RunAll.write(Line)
                            
                # loop to check if all simulations are complete
                RunAll.write('\n:: loop to wait till simulations are done\n')
                RunAll.write('SETLOCAL\n')
                RunAll.write(':WAIT\n')
                RunAll.write('TIMEOUT 10 > NUL\n')
                RunAll.write('SET /a COUNT=0\n')
                RunAll.write('FOR /F %%N IN (\'dir/s/b/a-d \"{:s}\\*.tmp'.format(TempDir) + \
                        '\"^|findstr /ric:"\\\\*\\.tmp$"^|find /c /v \"\"\') do set COUNT=%%N\n')
                RunAll.write('IF %COUNT% NEQ {:d} GOTO :WAIT\n'.format(iEnd-iStart))
                
                # Delete all temp files
                RunAll.write('\n:: delete all temp files\n')
                for iID in range(iStart,iEnd):
                    RunAll.write('DEL \"{:s}\\Done{:s}.tmp\"\n'.format(TempDir,FileIDs[iID]))
                    
                RunAll.write('\n:: run next bat file or state complete\n')
                if GrpID < NumGrps*len(TurbNames) - 1:
                    RunAll.write('@ECHO Calling batch file {:s}_{:s}_Group{:d}.bat\n'.format(RunName,ExeID,GrpID+1))
                    RunAll.write('{:s}\\{:s}_{:s}_Group{:d}.bat\n'.format(BaseBatDir,RunName,ExeID,GrpID+1))
                else:
                    RunAll.write('@ECHO All runs for {:s} complete\n'.format(ExeID))
                    
            GrpID += 1


