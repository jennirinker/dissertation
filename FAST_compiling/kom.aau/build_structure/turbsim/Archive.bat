@REM===================================================================================
@REM  BONNIE: CHECK DISCLAIMER TO THE ARCHIVE, ADDED PER LEGAL'S "REQUEST" VIA MARSHALL
@REM===================================================================================

@ECHO ON

@IF NOT "%1"==""  GOTO SetVars

@ECHO  
@ECHO  The syntax for creating an archive is "Archive <version> <AllFiles>"
@ECHO.
@ECHO  Example:  "archive 121"

@GOTO Done

@REM=================================================================================================================================
:SetVars
@SET ARCHNAME=TurbSim_v%1
@SET XTRANAME=TurbSim_v%1_ExtraFiles
@SET PROGNAME=TurbSim
@SET FILELIST=ArcFiles.txt
@SET XTRALIST=ArcXtraFiles.txt
@SET WINZIP="C:\Program Files (x86)\WinZip\WZZip"
@SET WINZIPSE="C:\Program Files (x86)\WinZip Self-Extractor\WZIPSE22\wzipse32.exe"

@IF EXIST ARCHTMP.zip DEL ARCHTMP.zip
@IF EXIST %ARCHNAME%.exe DEL %ARCHNAME%.exe

@SET ARCHPATH=Archive
@IF NOT "%3"=="" @SET ARCHPATH=Archive\AlphaVersions

@REM================================================================================================================================
@REM Create an archive of TurbSim in the format that is distributed on the web.                                             
@REM================================================================================================================================
:DoIt
@%WINZIP% -a -o -P ARCHTMP @%FILELIST%
@ECHO Creating self-extracting archive...
@%WINZIPSE% ARCHTMP.zip -d. -y -win32 -le -overwrite -st"Unzipping %PROGNAME%" -m Disclaimer.txt
@ECHO The self-extracting archive has been created.
@COPY ARCHTMP.exe %ARCHPATH%\%ARCHNAME%.exe
@DEL ARCHTMP.zip, ARCHTMP.exe


@IF NOT "%2"=="AllFiles" GOTO Done
@REM================================================================================================================================
@REM Create an archive of the "Extra" TurbSim documents that don't need to be released.  "AllFiles" 
@REM================================================================================================================================
@IF EXIST %XTRANAME%.zip DEL %XTRANAME%.zip  

@%WINZIP% -a -o -P ARCHTMP @%XTRALIST%                                      
@ECHO Creating self-extracting archive...                                   
@%WINZIPSE% ARCHTMP.zip -d. -y -win32 -le -overwrite -st"Unzipping %PROGNAME Extra Files" -m Disclaimer.txt
@ECHO The self-extracting archive has been created.                         
@COPY ARCHTMP.exe %ARCHPATH%\%XTRANAME%.exe                                 
@DEL ARCHTMP.zip, ARCHTMP.exe                                               


@REM================================================================================================================================
:Done
@REM================================================================================================================================
@SET ARCHNAME=
@SET XTRANAME=
@SET PROGNAME=
@SET FILELIST=
@SET XTRALIST=
@SET WINZIP=
@SET WINZIPSE=
