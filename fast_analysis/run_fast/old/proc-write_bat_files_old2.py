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
URefs  = [7.0, 9.0, 13.0]
Is     = [0.3]
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
                    
# -----------------------------------------------------------------------------

print('\nDeleting base directories...\n')

# clear BaseBatDir
if os.path.exists(BaseBatDir): shutil.rmtree(BaseBatDir)
os.mkdir(BaseBatDir)
if os.path.exists(BaseFastDir): shutil.rmtree(BaseFastDir)
os.mkdir(BaseFastDir)
if os.path.exists(BaseWindDir): shutil.rmtree(BaseWindDir)
os.mkdir(BaseWindDir)
    
# paths to group and individual bat templates
IndFullRunBatTmplPath = os.path.join(TmplDir,'Template_IndFullRunBat.bat')
TurbSimInpBatTmplPath = os.path.join(TmplDir,'Template_TurbSimInpBat.bat')
TurbSimBatTmplPath    = os.path.join(TmplDir,'Template_TurbSimBat.bat')
FastInpBatTmplPath    = os.path.join(TmplDir,'Template_FastInpBat.bat')
FastBatTmplPath       = os.path.join(TmplDir,'Template_FastBat.bat')

# loop through turbines
for TurbName in TurbNames:
    
    print('Writing bat files for turbine \"{:s}\"'.format(TurbName))
        
    # define folder for wind and FAST files
    BatDir  = os.path.join(BaseBatDir,TurbName)
    WindDir = os.path.join(BaseWindDir,TurbName)
    FastDir = os.path.join(BaseFastDir,TurbName)
    ModlDir = os.path.join(BaseModlDir,TurbName)

    # create directories
    os.mkdir(BatDir)
    os.mkdir(WindDir)
    os.mkdir(FastDir)
    MsgDir = os.path.join(BatDir,'Messages')
    os.mkdir(MsgDir)
    
    # copy pitch file to FastDir
    PitchStart = os.path.join(ModlDir,'pitch.ipt')
    PitchEnd   = os.path.join(FastDir,'pitch.ipt')
    shutil.copyfile(PitchStart,PitchEnd)
    
    # open RunAll file
    MonssonRunAllPath = os.path.join(BatDir,'{:s}_MonsoonRunAll.bat'.format(RunName))
    RunAllPath = os.path.join(BatDir,'{:s}_RunAll.bat'.format(RunName))
    with open(MonssonRunAllPath,'w') as MonsoonRunAll:
      with open(RunAllPath,'w') as RunAll:
          
        # unmap used drives
        MonsoonRunAll.write('net use A: /delete\n')
        RunAll.write('net use A: /delete\n')
        MonsoonRunAll.write('net use S: /delete\n')
        RunAll.write('net use S: /delete\n')
        
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
                
                # define paths for bat files and messages
                IndFullRunBatName = '{:s}_IndFullRun.bat'.format(FileID)
                IndFullRunBatPath = os.path.join(BatDir,IndFullRunBatName)
                IndFullRunBatMsgName = os.path.splitext( \
                                            IndFullRunBatName)[0] + '.mes'
                IndFullRunBatMsgPath = os.path.join(MsgDir,IndFullRunBatMsgName)
                
                TurbSimInpBatName = '{:s}_TurbSimInp.bat'.format(FileID)
                TurbSimInpBatPath = os.path.join(BatDir,TurbSimInpBatName)
                TurbSimInpBatMsgName = os.path.splitext( \
                                            TurbSimInpBatName)[0] + '.mes'
                TurbSimInpBatMsgPath = os.path.join(MsgDir,TurbSimInpBatMsgName)
                
                TurbSimBatName = '{:s}_TurbSim.bat'.format(FileID)
                TurbSimBatPath = os.path.join(BatDir,TurbSimBatName)
                TurbSimBatMsgName = os.path.splitext( \
                                            TurbSimBatName)[0] + '.mes'
                TurbSimBatMsgPath = os.path.join(MsgDir,TurbSimBatMsgName)
                
                FastInpBatName = '{:s}_FastInp.bat'.format(FileID)
                FastInpBatPath = os.path.join(BatDir,FastInpBatName)
                FastInpBatMsgName = os.path.splitext( \
                                            FastInpBatName)[0] + '.mes'
                FastInpBatMsgPath = os.path.join(MsgDir,FastInpBatMsgName)
                
                FastBatName = '{:s}_Fast.bat'.format(FileID)
                FastBatPath = os.path.join(BatDir,FastBatName)
                FastBatMsgName = os.path.splitext( \
                                            FastBatName)[0] + '.mes'
                TurbSimInpBatMsgPath = os.path.join(MsgDir,FastBatMsgName)
                
                # define paths for TurbSim and Fast files
                TurbSimInpName = '{:s}_{:s}.inp'.format(TurbName,FileID)
                TurbSimName    = '{:s}_{:s}.bts'.format(TurbName,FileID)
                FastInpName    = '{:s}_{:s}.fst'.format(TurbName,FileID)
                FastOutName    = '{:s}_{:s}.out'.format(TurbName,FileID)
                TurbSimInpPath = os.path.join(WindDir,TurbSimInpName)
                TurbSimOutPath = os.path.join(WindDir,TurbSimName)
                FastInpPath    = os.path.join(FastDir,FastInpName)
                FastOutPath    = os.path.join(FastDir,FastOutName)
        
                # create time stamp comment
                TimeStr   = datetime.datetime.now().strftime(DateFmt)
                IndFullRunBatCmnt = '::   created {:s}'.format(TimeStr) + \
                                    ' by J. Rinker Python code'
        
                # open the full run bat file, loop through lines
                with open(IndFullRunBatTmplPath,'r') as Tmpl:
                    with open(IndFullRunBatPath,'w') as Bat:
                        for Line in Tmpl:
                            
                            # if line has a writeable field
                            if ('{:s' in Line):
                                Field   = Line.split('{:s<')[1].split('>}')[0]
                                NewLine = ''.join(Line.split( \
                                            '<{:s}>'.format(Field))).format(eval(Field))
                                
                            else:
                                NewLine = Line
                                
                            Bat.write(NewLine)
                                
                # open the turbsim input bat file, loop through lines
                TurbSimInpCmnd = '{:s} {:s} {:s} {:s}'.format(TurbName,
                                                        RunName,FileID,ParmStr)
                with open(TurbSimInpBatTmplPath,'r') as Tmpl:
                    with open(TurbSimInpBatPath,'w') as Bat:
                        for Line in Tmpl:
                            
                            # if line has a writeable field
                            if ('{:s' in Line):
                                Field   = Line.split('{:s<')[1].split('>}')[0]
                                NewLine = ''.join(Line.split( \
                                            '<{:s}>'.format(Field))).format(eval(Field))
                                
                            else:
                                NewLine = Line
                                
                            Bat.write(NewLine)
                                
                # open the turbsim bat file, loop through lines
                with open(TurbSimBatTmplPath,'r') as Tmpl:
                    with open(TurbSimBatPath,'w') as Bat:
                        for Line in Tmpl:
                            
                            # if line has a writeable field
                            if ('{:s' in Line):
                                Field   = Line.split('{:s<')[1].split('>}')[0]
                                NewLine = ''.join(Line.split( \
                                            '<{:s}>'.format(Field))).format(eval(Field))
                                
                            else:
                                NewLine = Line
                                
                            Bat.write(NewLine)
                                
                # open the fast input bat file, loop through lines
                FastInpCmnd = '{:s} {:s} {:s}'.format(TurbName,RunName,
                                                      TurbSimOutPath)
                with open(FastInpBatTmplPath,'r') as Tmpl:
                    with open(FastInpBatPath,'w') as Bat:
                        for Line in Tmpl:
                            
                            # if line has a writeable field
                            if ('{:s' in Line):
                                Field   = Line.split('{:s<')[1].split('>}')[0]
                                NewLine = ''.join(Line.split( \
                                            '<{:s}>'.format(Field))).format(eval(Field))
                                
                            else:
                                NewLine = Line
                                
                            Bat.write(NewLine)
                                
                # open the fast bat file, loop through lines
                with open(FastBatTmplPath,'r') as Tmpl:
                    with open(FastBatPath,'w') as Bat:
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
                        ' /jobname:{:s} \"{:s}"\n'.format( \
                                IndFullRunBatName,IndFullRunBatPath)
                MonsoonRunAll.write(Line)
                Line = 'start /b {:s}\n'.format(IndFullRunBatPath)
                RunAll.write(Line)



