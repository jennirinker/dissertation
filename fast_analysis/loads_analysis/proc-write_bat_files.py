"""
Write bat file for monsoon simulation
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\public\\nwtc_python_tools'
if (libpath not in sys.path): sys.path.append(libpath)

import os
import datetime
import shutil
import random

# define turbine name, turbine dictionary
TurbNames = ['WP0.75A08V00','WP1.5A08V03','WP3.0A02V02','WP5.0A04V00']
#TurbNames = ['WP0.75A08V00']


# run specifics
DateFmt = '%a %b %d %H:%M:%S %Y'            # date format
RunName = 'SmallRun'                        # run name
URefs  = [7.0, 9.0, 10.0, 11.0, 13.0, 19.0]
Is     = [0.1,0.3,0.5]
logLs  = [2.0]
rhos   = [0.4]
n_dups = 10

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
#TurbSimInpExePath = os.path.join(ExeDir,'TurbSimInp.py')
#TurbSimExePath = os.path.join(ExeDir,'TurbSim.exe')
#FastInpExePath = os.path.join(ExeDir,'FastInp.py')
#FastExePath    = os.path.join(ExeDir,'FAST.exe')
TurbSimInpExePath = 'S:\\\\TurbSimInp.py'
TurbSimExePath = 'S:\\\\TurbSim.exe'
FastInpExePath = 'S:\\\\FastInp.py'
FastExePath    = 'S:\\\\FAST.exe'
ExeNames = ['TurbSimInp.py','TurbSim.exe','FastInp.py','FAST.exe']
ExeIDs   = ['TSInp','RunTS','FSInp','RunFS']
                    
# -----------------------------------------------------------------------------

print('\nDeleting base directories...\n')

# clear BaseBatDir
if os.path.exists(BaseBatDir): shutil.rmtree(BaseBatDir)
os.mkdir(BaseBatDir)
if os.path.exists(BaseFastDir): shutil.rmtree(BaseFastDir)
os.mkdir(BaseFastDir)
if os.path.exists(BaseWindDir): shutil.rmtree(BaseWindDir)
os.mkdir(BaseWindDir)
    
# paths to bat template
BatTmplPath = os.path.join(TmplDir,'Template.bat')

# loop through TSInp, RunTS, FSInp, and RunFS
for iExe in range(len(ExeNames)):
    
    ExeName = ExeNames[iExe]
    ExeID   = ExeIDs[iExe]
    ExePath = os.path.join(ExeDir,ExeName)
    TimeStr   = datetime.datetime.now().strftime(DateFmt)
    BatCmnt = '::   bat file to run {:s} '.format(ExeName) + \
                 'created {:s} by J. Rinker Python code'.format(TimeStr)

    print('Writing bat files for \"{:s}\"'.format(ExeID))
    
    # open RunAll file
    RunAllPath = os.path.join(BaseBatDir,
                              '{:s}_{:s}_RunAll.bat'.format(RunName,ExeID))
                              
    with open(RunAllPath,'w') as RunAll:
    
        # unmap used drives
        RunAll.write('net use A: /delete\n')
        RunAll.write('net use S: /delete\n')
        
        # loop through turbines
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
    
            # iterate over all wind parameters
            for iU,iI,iL,iR in [(a,b,c,d) for a in range(len(URefs)) \
                                        for b in range(len(Is)) \
                                        for c in range(len(logLs)) \
                                        for d in range(len(rhos))]:
                URef,I,L,rho = URefs[iU],Is[iI],10**logLs[iL],rhos[iR]
                                            
                # loop through duplicates               
                for iS in range(n_dups):
    
                    # generate random seeds
                    R1 = random.randint(-2147483648,2147483647)
                    R2 = random.randint(-2147483648,2147483647)
    
                    # create parameter string
                    ParmStr = '\"{:.2f} {:.3f} {:.1f}'.format(\
                                URef,I*URef,L)  + \
                        ' {:.2f} {:d} {:d}\"'.format(\
                            rho,R1,R2)
    
                    # create file ID
                    FileID = format(iU,'x')+format(iI,'x')+ \
                                format(iL,'x')+format(iR,'x')+ \
                                format(iS,'x') 
                
                    # bat and message file names and paths
                    BatName    = '{:s}_{:s}_{:s}.bat'.format(TurbName,FileID,ExeID)
                    BatMsgName = '{:s}_{:s}_{:s}.mes'.format(TurbName,FileID,ExeID)
                    BatPath    = os.path.join(BatDir,BatName)
                    
                    # paths to TurbSim and FAST input and output files
                    TurbSimInpName = '{:s}_{:s}.inp'.format(TurbName,FileID)
                    TurbSimOutName = '{:s}_{:s}.bts'.format(TurbName,FileID)
                    FastInpName    = '{:s}_{:s}.fst'.format(TurbName,FileID)
                    TurbSimInpPath = os.path.join(WindDir,TurbSimInpName)
                    TurbSimOutPath = os.path.join(WindDir,TurbSimOutName)
                    FastInpPath    = os.path.join(FastDir,FastInpName)
        
                    # exectuable command
                    if (iExe == 0):
                        ExeCmnd = 'python %EXE% ' + \
                                    '{:s} {:s} {:s} {:s}'.format(TurbName,
                                                        RunName,FileID,ParmStr)
                    elif (iExe == 1):
                        ExeCmnd = '%EXE% {:s}'.format(TurbSimInpPath)
                    elif (iExe == 2):
                        ExeCmnd = 'python %EXE% {:s} {:s} {:s}'.format( \
                                    TurbName,RunName,TurbSimOutPath)
                    elif (iExe == 3):
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
                    if (iExe % 2):
                        Line = 'start /b job submit /scheduler:monsoon' + \
                            ' /jobname:{:s} \"{:s}"\n'.format(BatName,BatPath)
                    else:
                        Line = 'start /b {:s}\n'.format(BatPath)
                    RunAll.write(Line)
                    
        RunAll.write('@ECHO All runs for {:s} complete\n'.format(ExeID))


