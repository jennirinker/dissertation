@ECHO OFF

set LOC=Y:\Wind\WindWeb\designcodes\preprocessors\turbsim

@IF /I "%1"=="-NOFILES" GOTO UpdateWeb
@IF NOT "%1"=="" GOTO CopyFiles

@ECHO
@ECHO  The syntax for updating files on the TurbSim page is "update_web <version>"
@ECHO.
@ECHO  Example:  "update_web 1.60.00"
@ECHO.
@ECHO  The syntax for updating the web page without copying the other files is "update_web -nofiles"
@ECHO.

@GOTO Done

:CopyFiles
copy ..\ChangeLog.txt                  %LOC%\ChangeLog.txt
copy ..\Archive\TurbSim_v%1.exe        %LOC%\TurbSim_v%1.exe
copy ..\TurbSim.pdf                    %LOC%\TurbSim.pdf

rem THESE FILES SHOULDN'T HAVE TO BE UPDATED AGAIN
REM copy ..\TurbSimOverview.pdf            %LOC%\TurbSimOverview.pdf
REM copy M:\coh_events\Archive\TSM_structures.exe        %LOC%\TSM_structures.exe
REM copy M:\coh_events\Archive\TSM_DNS_structures.exe    %LOC%\TSM_DNS_structures.exe

:UpdateWeb
copy    CreatePage.pl                  %LOC%\CreatePage.pl
perl %LOC%\CreatePage.pl


:Done
set LOC=

