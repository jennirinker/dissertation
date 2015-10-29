REM Turn echo off to prevent "REM" lines from echoing to command line
@echo off

REM
REM Windows 7 batch file to automate unzipping of files for compiling
REM FAST, AeroDyn, MCrunch, and TurbSim.
REM
REM Written by Jenni Rinker, 19-October-2014.
REM
REM Modified from http://kom.aau.dk/~anb/fast/fastinstall.html
REM
REM REQUIREMENTS:
REM  - Programs (installed and added to system PATH):
REM	Cygwin
REM	7-zip
REM  - Files (from URL above, in same directory as .bat):
REM	build_structure.tar.gz
REM	AD_v13.00.02a-bjj.tar.gz
REM	FAST_v7.02.00d-bjj.tar.gz
REM	NWTC_Lib_v1.07.00b-mlb.tar.gz
REM	Crunch_v3.02.00c-mlb.tar.gz
REM	InflowWind_v1.02.00b-adp.exe
REM	MBC_v1.00.00a-gsb.exe
REM	MCrunch_v1.00.00ab-gjh.exe
REM	TurbSim_v1.06.00.exe

REM Choose directory where files are located
set "dir_path=C:\Users\jrinker\Dropbox\research_code\kom.aau"

REM Print chosen directory to screen
echo Directory path: %dir_path%

REM Change to that directory
cd %dir_path%

REM Set directories and option for extracting with 7-zip
set "inf_opt=-o%dir_path%\build_structure\inflow"
set "mbc_opt=-o%dir_path%\build_structure\mbc"
set "mcr_opt=-o%dir_path%\build_structure\mcrunch"
set "tur_opt=-o%dir_path%\build_structure\turbsim"

REM Use Cygwin to unzip tar.gz files to their respective subdirectories
REM  NOTE: relative pathnames must have forwardslash for Cygwin to unzip
tar -zxvf build_structure.tar.gz
tar -zxvf AD_v13.00.02a-bjj.tar.gz -C build_structure/aero/
tar -zxvf FAST_v7.02.00d-bjj.tar.gz -C build_structure/fast/
tar -zxvf NWTC_Lib_v1.07.00b-mlb.tar.gz -C build_structure/nwtc/
tar -zxvf Crunch_v3.02.00c-mlb.tar.gz -C build_structure/crunch/

REM Use 7-zip to extract files to specific folders
7z x InflowWind_v1.02.00b-adp.exe %inf_opt%
7z x MBC_v1.00.00a-gsb.exe %mbc_opt%
7z x MCrunch_v1.00.00ab-gjh.exe %mcr_opt%
7z x TurbSim_v1.06.00.exe %tur_opt%
