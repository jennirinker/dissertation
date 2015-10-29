This archive is designed to test the InflowWind module, scheduled for release in March 2010.

To compile, you will need to add the following files to your project
(1) Files from NWTC_Library v1.04.01:
      SingPrec.f90 ( or DoubPrec.f90 )
      SysIVF.f90 ( or another Sys*.f90 file, depending on the system you're compiling on)
      NWTC_IO.f90
      NWTC_Num.f90
      NWTC_Aero.f90
      NWTC_Library.f90

(2) All the files from the InflowWind\Source folder:
      SharedInflowDefs.f90
      CTWind.f90
      FDWind.f90
      HAWCWind.f90
      FFWind.f90
      HHWind.f90
      UserWind.f90
      InflowWindMod.f90

(3) The Fortran test program:
      InflowWind.f90   
      
The files should be compiled in the order specified above.