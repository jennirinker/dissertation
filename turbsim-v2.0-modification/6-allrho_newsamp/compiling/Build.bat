@ECHO OFF

REM ----------------------------------------------------------------------------
REM                   set compiler internal variables
REM ----------------------------------------------------------------------------
REM    You can run this bat file from the IVF compiler's command prompt (and not
REM    do anything in this section). If you choose not to run from the IVF command
REM    prompt, you must call the compiler's script to set internal variables.
REM    TIP: Right click on the IVF Compiler's Command Prompt shortcut, click
REM    properties, and copy the target (without cmd.exe and/or its switches) here:

rem IF "%INTEL_SHARED%"=="" CALL "C:\Program Files\Intel\Compiler\Fortran\10.1.024\IA32\Bin\IFORTVARS.bat"
call "C:\Program Files (x86)\Intel\Composer XE 2011 SP1\bin\ipsxe-comp-vars.bat" ia32 vs2010
rem CALL "C:\Program Files (x86)\Intel\Composer XE 2011 SP1\bin\ipsxe-comp-vars.bat" intel64 vs2010



REM ----------------------------------------------------------------------------
REM -------------------- LOCAL VARIABLES ---------------------------------------
REM ----------------------------------------------------------------------------

@SET ROOT_NAME=TurbSim
@SET EXEname=%ROOT_NAME%.exe

@SET COPTS=/threads /O2 /inline:speed /traceback /real_size:32 /assume:byterecl /nostand
@SET LOPTS=/link /stack:64000000


@REM ----------------------------------------------------------------------------
@REM ------------------------- LOCAL PATHS --------------------------------------
@REM ----------------------------------------------------------------------------
@REM -- USERS WILL NEED TO EDIT THESE PATHS TO POINT TO FOLDERS ON THEIR LOCAL --
@REM -- MACHINES.  NOTE: do not use quotation marks around the path names!!!! ---
@REM ----------------------------------------------------------------------------
@REM NWTC_Lib_Loc is the location of the NWTC subroutine library files
@REM TurbSim_Loc  is the location of the TurbSim (main) source files
@REM ----------------------------------------------------------------------------

@SET NWTC_Lib_Loc=C:\Users\bjonkman\Documents\DATA\DesignCodes\miscellaneous\nwtc_subs\SVNdirectory\trunk_V1\source
@SET TurbSim_Loc=C:\Users\bjonkman\Documents\DATA\DesignCodes\preprocessors\TurbSim\SVNdirectory\branches\APIABSWnd\Source



@REM ----------------------------------------------------------------------------
@REM -------------------- LIST OF ALL SOURCE FILES ------------------------------
@REM ----------------------------------------------------------------------------
:SourceFiles


@REM BJJ: Uses v1.07.02, 18-Jan-2012

@SET NWTC_Files="%NWTC_Lib_Loc%\SingPrec.f90" ^
                "%NWTC_Lib_Loc%\SysIVF.f90" ^
                "%NWTC_Lib_Loc%\NWTC_IO.f90" ^
                "%NWTC_Lib_Loc%\NWTC_Num.f90" ^
                "%NWTC_Lib_Loc%\NWTC_Aero.f90" ^
                "%NWTC_Lib_Loc%\ModMesh.f90" ^
                "%NWTC_Lib_Loc%\NWTC_Library.f90"

@SET TS_Lib_Files="%TurbSim_Loc%\RNG\RANLUX.f90" ^
                  "%TurbSim_Loc%\NetLib\FFTPACK\fftpack.f" ^
                  "%TurbSim_Loc%\NetLib\LAPACK\spptrf.f" ^
                  "%TurbSim_Loc%\NetLib\LAPACK\BLAS\BLAS2.f" ^
                  "%TurbSim_Loc%\NetLib\LAPACK\BLAS\sdot.f" ^
                  "%TurbSim_Loc%\NetLib\LAPACK\BLAS\sscal.f" ^
                  "%TurbSim_Loc%\NetLib\LAPACK\BLAS\xerbla.f"

@SET TS_Files="%TurbSim_Loc%\FFTMod.f90" ^
              "%TurbSim_Loc%\Modules.f90" ^
              "%TurbSim_Loc%\Root_Searching.f90" ^
              "%TurbSim_Loc%\BlankModVKM.f90" ^
              "%TurbSim_Loc%\TSsubs.f90" ^
              "%TurbSim_Loc%\SetVersion.f90" ^
              "%TurbSim_Loc%\TurbSim.f90"


@REM ----------------------------------------------------------------------------
@REM ---------------- COMPILE WITH INTEL VISUAL FORTRAN -------------------------
@REM ----------------------------------------------------------------------------
:ivf

ECHO.
ECHO Compiling TurbSim and NWTC Library routines to create %EXEname%:
ifort %COPTS% %NWTC_Files% %TS_Lib_Files% %TS_Files% %LOPTS% /out:%EXEname%


REM ----------------------------------------------------------------------------
REM ------------------------- CLEAR MEMORY -------------------------------------
REM ------------- and delete all .mod and .obj files ---------------------------
REM ----------------------------------------------------------------------------
:end
ECHO 

@DEL *.mod
@DEL *.obj

@SET ROOT_NAME=
@SET EXEname=

@SET COPTS=
@SET LOPTS=

@SET NWTC_Lib_Loc=
@SET TurbSim_Loc=


@SET NWTC_Files=
@SET TS_Lib_Files=
@SET TS_Files=





