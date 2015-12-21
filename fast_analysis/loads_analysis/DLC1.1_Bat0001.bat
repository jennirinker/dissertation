:: Remote Execution Batch File: \\wind-dfs1\Public\Projects\Projects T-Z\Tip Speed\LoadsAnalysis\Baseline--original_blades\LoadsAnalysis\\DLC1.1\DLC1.1_Bat0001.bat
::   created Thu Oct 31 15:17:44 2013 by ..\..\Software\RunIEC\RunIEC.pl

::@ECHO Off

@NET USE A:    "\\wind-dfs1\Public\Projects\Projects T-Z\Tip Speed\LoadsAnalysis\Baseline--original_blades\LoadsAnalysis"    /persistent:no
@SET MESSFILE="A:\\DLC1.1\Messages\DLC1.1_Bat0001.mes"
@ECHO Mounting file servers.                                         >  %MESSFILE%
@NET USE N:       "\\wind-dfs1\Public\Projects\Projects T-Z\Tip Speed\LoadsAnalysis\Baseline--original_blades\MasterModel"       /persistent:no >> %MESSFILE%
@NET USE S:    "\\wind-dfs1\Public\Projects\Projects T-Z\Tip Speed\LoadsAnalysis\Software"    /persistent:no >> %MESSFILE%
@NET USE P: "\\wind-dfs1\Public\Projects\Projects T-Z\Tip Speed\LoadsAnalysis\Environment" /persistent:no >> %MESSFILE%
@NET USE R:    "\\wind-dfs1\Public\Projects\Projects T-Z\Tip Speed\LoadsAnalysis\Baseline--original_blades\MasterModel\AeroData"    /persistent:no >> %MESSFILE%
@NET USE W:       "\\wind-dfs1\Public\Projects\Projects T-Z\Tip Speed\LoadsAnalysis\Baseline--original_blades\MasterModel"       /persistent:no >> %MESSFILE%

@rem #mlb NET USE

@SET HostName=%COMPUTERNAME%
@SET FAST=S:\FAST\FAST7.02_AD13.00.02_MoreStrain_DLL.exe
@IF EXIST %FAST% GOTO FASTOK
  @ECHO Fast executable '\\wind-dfs1\Public\Projects\Projects T-Z\Tip Speed\LoadsAnalysis\Software\FAST\FAST7.02_AD13.00.02_MoreStrain_DLL.exe' does not exist!
  @ECHO RunIEC job is terminating
  @ECHO Fast executable '\\wind-dfs1\Public\Projects\Projects T-Z\Tip Speed\LoadsAnalysis\Software\FAST\FAST7.02_AD13.00.02_MoreStrain_DLL.exe' does not exist! >> %MESSFILE%
  @ECHO RunIEC job is terminating >> %MESSFILE%
  @EXIT
:FASTOK

:: Save current drive so we go back to the right place
FOR %%A in ("%CD%") DO SET MyDrive=%%~dA

@ECHO This job is running on: %HostName% at %Date% %Time%   >> %MESSFILE%
@ECHO   Batch file \\wind-dfs1\Public\Projects\Projects T-Z\Tip Speed\LoadsAnalysis\Baseline--original_blades\LoadsAnalysis\\DLC1.1\DLC1.1_Bat0001.bat                            >> %MESSFILE%
@ECHO   FAST: %FAST%                                        >> %MESSFILE%

:: Move to A:

@NET USE             >> %MESSFILE%
@A:   >> %MESSFILE%
FOR %%A in ("%CD%") DO SET RunDrive=%%~dA
@rem #mlb ECHO Running from: %RunDrive% (\\wind-dfs1\Public\Projects\Projects T-Z\Tip Speed\LoadsAnalysis\Baseline--original_blades\LoadsAnalysis)
@rem #mlb ECHO Running from: %RunDrive% (\\wind-dfs1\Public\Projects\Projects T-Z\Tip Speed\LoadsAnalysis\Baseline--original_blades\LoadsAnalysis) >> %MESSFILE%

:: =================================== Simulations ======================

:: --- Simulation 1 ---
@ECHO  Running SIM 1 DLC1.1\DLC1.1_0001_Land_04.0V0_S01.bat
@ECHO  Running SIM 1 DLC1.1\DLC1.1_0001_Land_04.0V0_S01.bat >> %MESSFILE%
@CALL "A:\\DLC1.1\DLC1.1_0001_Land_04.0V0_S01.bat"
@ECHO    SIM 1  ERRORLEVEL %ERRORLEVEL% >> %MESSFILE%



::%MyDrive% 2>&1
