"""
Write bat file for monsoon simulation: just TSRun and FSRun
"""
import sys
libpath = 'C:\\Users\\jrinker\\Documents\\GitHub\\dissertation'
if (libpath not in sys.path): sys.path.append(libpath)

import JR_Library.main as jr
import os
import datetime
import shutil
import numpy as np

# define turbine name, turbine dictionary
#TurbNames = ['WP0.75A08V00','WP1.5A08V03','WP3.0A02V02','WP5.0A04V00']
TurbNames = ['WP5.0A04V00']


# run specifics
DateFmt    = '%a %b %d %H:%M:%S %Y'            # date format
SimsPerBat = 100
RunName = 'Fine'                        # run name
RunParms = jr.RunName2WindParms(RunName)
URefs,Is,Ls,rhos,n_dups  = RunParms['URefs'],RunParms['Is'],RunParms['Ls'], \
                            RunParms['rhos'],RunParms['n_dups']

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
if os.path.exists(BaseFastDir): shutil.rmtree(BaseFastDir)
os.mkdir(BaseFastDir)
    
# paths to bat template
BatTmplPath = os.path.join(TmplDir,'Template.bat')

# number of bat files
NumSims = len(URefs)*len(Is)*len(Ls)*len(rhos)*n_dups
NumBatFiles = int(np.ceil(NumSims/float(SimsPerBat)))

# create list of file ids
FileIDs = []
for iU,iI,iL,iR in [(a,b,c,d) for a in range(len(URefs)) \
                            for b in range(len(Is)) \
                            for c in range(len(Ls)) \
                            for d in range(len(rhos))]:
                                
    # loop through duplicates               
    for iS in range(n_dups):

        # create file ID
        FileIDs.append(format(iU,'x')+format(iI,'x')+ \
                    format(iL,'x')+format(iR,'x')+ \
                    format(iS,'x'))

# loop through executable names
for iExe in range(len(ExeNames)):
    
    ExeName = ExeNames[iExe]
    ExeID   = ExeIDs[iExe]
    ExePath = os.path.join(ExeDir,ExeName)
    TimeStr   = datetime.datetime.now().strftime(DateFmt)
    BatCmnt = '::   bat file to run {:s} '.format(ExeName) + \
                 'created {:s} by J. Rinker Python code'.format(TimeStr)

    print('Writing bat files for \"{:s}\"'.format(ExeID))
    
    # loop through RunAll bat files
    for iBat in range(NumBatFiles):
    
        # open RunAll file
        RunAllPath = os.path.join(BaseBatDir,
                                  '{:s}_{:s}_RunAll{:d}.bat'.format(\
                                                      RunName,ExeID,iBat))
        iSim = 0
                                  
        with open(RunAllPath,'w') as RunAll:
        
            # unmap used drives
            RunAll.write('net use A: /delete\n')
            RunAll.write('net use S: /delete\n')
            RunAll.write('@ECHO Off\n')
            
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
                while iSim < SimsPerBat:
                    
                    FileID = FileIDs[iSim+iBat*SimsPerBat]
                    
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
                    
                    iSim += 1
                        


