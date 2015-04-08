@ECHO OFF

@IF NOT "%1"==""  GOTO Copy

@ECHO  
@ECHO  The syntax for updating files on the TurbSim alpha page is "update_web <version>"
@ECHO.
@ECHO  Example:  "update_web 120"

@GOTO Done

:Copy

set LOC=Y:\Wind\WindWeb\designcodes\preprocessors\turbsim\alpha
set Alpha_Dir=D:\DATA\DesignCodes\preprocessors\TurbSim\SVNdirectory\branches\PeriodicOutput


copy %Alpha_Dir%\ChangeLog.txt                  %LOC%\ChangeLog.txt
copy CreatePage_alpha.pl                        %LOC%\CreatePage.pl
copy %Alpha_Dir%\Archive\TurbSim_v%1.exe        %LOC%\TurbSim_v%1.exe
REM copy ..\TurbSim.pdf                        %LOC%\TurbSim.pdf
REM copy ..\TurbSimOverview.pdf            %LOC%\TurbSimOverview.pdf

REM copy M:\coh_events\Archive\TSM_structures.exe        %LOC%\TSM_structures.exe
REM copy M:\coh_events\Archive\TSM_DNS_structures.exe    %LOC%\TSM_DNS_structures.exe

perl %LOC%\CreatePage.pl

set LOC=

:Done
