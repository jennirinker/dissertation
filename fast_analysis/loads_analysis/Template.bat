
:: ====================== {:s} ========================

@SET SMSSFILE="{:s}.mes"
@ECHO  ========= {:s} ========= >  %SMSSFILE%
@ECHO  Running FAST for case {:s}  >> %SMSSFILE%

%FAST% {:s}                   >> %SMSSFILE%
   IF ERRORLEVEL 1 GOTO :ERROR1
@ECHO Case {:s} successfully ran.      >> %SMSSFILE%
@ECHO  Output is in {:s} (under {:s}). >> %SMSSFILE%
GOTO :SIMS1

:ERROR1
@ECHO   ==( Case {:s} failed )==
@ECHO   ==( Case {:s} failed )== >> %SMSSFILE%

:SIMS1
@ECHO This job ran on: %HostName% finishing at %Date% %Time%  >> %SMSSFILE%
:: EXIT /b 0
