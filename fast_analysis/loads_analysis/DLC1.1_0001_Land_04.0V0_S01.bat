
:: =================================== Simulation 1 ======================

@SET SMSSFILE="A:\\DLC1.1\Messages\DLC1.1_0001_Land_04.0V0_S01.mes"
@ECHO  ========= Simulation 1 >  %SMSSFILE%
@ECHO  Running FAST for case \DLC1.1\DLC1.1_0001_Land_04.0V0_S01.  >> %SMSSFILE%

%FAST% \DLC1.1\DLC1.1_0001_Land_04.0V0_S01.fst                   >> %SMSSFILE%
   IF ERRORLEVEL 1 GOTO :ERROR1
@ECHO Case DLC1.1_0001_Land_04.0V0_S01 successfully ran.      >> %SMSSFILE%
@ECHO  Output is in \DLC1.1\DLC1.1_0001_Land_04.0V0_S01.out (under \\wind-dfs1\Public\Projects\Projects T-Z\Tip Speed\LoadsAnalysis\Baseline--original_blades\LoadsAnalysis). >> %SMSSFILE%
GOTO :SIMS1

:ERROR1
@ECHO   ==( Case DLC1.1_0001_Land_04.0V0_S01 failed )==
@ECHO   ==( Case DLC1.1_0001_Land_04.0V0_S01 failed )== >> %SMSSFILE%

:SIMS1
@ECHO This job ran on: %HostName% finishing at %Date% %Time%  >> %SMSSFILE%
:: EXIT /b 0
