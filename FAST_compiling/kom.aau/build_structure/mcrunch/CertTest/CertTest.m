% Run MCrunch for a series of sample settings files.
%
% Syntax is:  CertTest
%
% Example:
%     CertTest
%
% See also CompFatigue, FileComp, GenPDFs, GenPSDs, GenRFPlots, GenStats, GenTimePlots, GetRoot,
%          MCrunch, ReadManyFiles, ReadSettings, ReadVal


   % Delete the old comparison file.  Add a header.

DelFile( 'CertTest.comp' )

fid = fopen( 'CertTest.comp', 'wt' );

if ( fid < 0 )
   beep
   error( '  Could not open "CertTest.comp" for reading.' );
end

DateTime = clock;
Date     = date;

fprintf( fid, '\nThis certification comparison was generated on %s at %02d:%02d:%02d.\n', Date, uint8( DateTime(4:6) ) );

fclose( fid );


fprintf( '\n=======\nTest_01\n=======\n' );
MCrunch 'Test_01.mcru'
if ( ~exist( 'Test01', 'dir' ) )
   mkdir( 'Test01' )
end % if
copyfile( 'DLC2.3_1.rflo'   , 'Test01\DLC2.3_1.rflo'    )
copyfile( 'RootFxc1.sums'   , 'Test01\RootFxc1.sums'    )
copyfile( 'RootFyc1.sums'   , 'Test01\RootFyc1.sums'    )
copyfile( 'RootFzc1.sums'   , 'Test01\RootFzc1.sums'    )
copyfile( 'Test_01.dels'    , 'Test01\Test_01.dels'     )
copyfile( 'Test_01_DEL.html', 'Test01\Test_01_DEL.html' )
FileComp( 'Test01\DLC2.3_1.rflo'   , 'TestFiles\Test01\DLC2.3_1.rflo'   , 5, 'CertTest.comp' )
FileComp( 'Test01\RootFxc1.sums'   , 'TestFiles\Test01\RootFxc1.sums'   , 5, 'CertTest.comp' )
FileComp( 'Test01\RootFyc1.sums'   , 'TestFiles\Test01\RootFyc1.sums'   , 5, 'CertTest.comp' )
FileComp( 'Test01\RootFzc1.sums'   , 'TestFiles\Test01\RootFzc1.sums'   , 5, 'CertTest.comp' )
FileComp( 'Test01\Test_01.dels'    , 'TestFiles\Test01\Test_01.dels'    , 5, 'CertTest.comp' )
FileComp( 'Test01\Test_01_DEL.html', 'TestFiles\Test01\Test_01_DEL.html', 5, 'CertTest.comp' )

fprintf( '\n=======\nTest_02\n=======\n' );
MCrunch 'Test_02.mcru'
if ( ~exist( 'Test02', 'dir' ) )
   mkdir( 'Test02' )
end % if
copyfile( 'Test_02_DEL.html', 'Test02\Test_02_DEL.html' )
FileComp( 'Test02\Test_02_DEL.html', 'TestFiles\Test02\Test_02_DEL.html', 5, 'CertTest.comp' )

fprintf( '\n=======\nTest_03\n=======\n' );
MCrunch 'Test_03.mcru'
if ( ~exist( 'Test03', 'dir' ) )
   mkdir( 'Test03' )
end % if
copyfile( 'DLC2.3_1.extr'  , 'Test03\DLC2.3_1.extr'   )
copyfile( 'DLC2.3_2.extr'  , 'Test03\DLC2.3_2.extr'   )
copyfile( 'RootFxc1.sums'  , 'Test03\RootFxc1.sums'   )
copyfile( 'RootFyc1.sums'  , 'Test03\RootFyc1.sums'   )
copyfile( 'RootFzc1.sums'  , 'Test03\RootFzc1.sums'   )
copyfile( 'Test03_Agg.extr', 'Test03\Test03_Agg.extr' )
FileComp( 'Test03\DLC2.3_1.extr'  , 'TestFiles\Test03\DLC2.3_1.extr'  , 5, 'CertTest.comp' )
FileComp( 'Test03\DLC2.3_2.extr'  , 'TestFiles\Test03\DLC2.3_2.extr'  , 5, 'CertTest.comp' )
FileComp( 'Test03\RootFxc1.sums'  , 'TestFiles\Test03\RootFxc1.sums'  , 5, 'CertTest.comp' )
FileComp( 'Test03\RootFyc1.sums'  , 'TestFiles\Test03\RootFyc1.sums'  , 5, 'CertTest.comp' )
FileComp( 'Test03\RootFzc1.sums'  , 'TestFiles\Test03\RootFzc1.sums'  , 5, 'CertTest.comp' )
FileComp( 'Test03\Test03_Agg.extr', 'TestFiles\Test03\Test03_Agg.extr', 5, 'CertTest.comp' )

fprintf( '\n=======\nTest_04\n=======\n' );
MCrunch 'Test_04.mcru'
if ( ~exist( 'Test04', 'dir' ) )
   mkdir( 'Test04' )
end % if
copyfile( 'DLC2.3_1.rflo', 'Test04\DLC2.3_1.rflo' )
copyfile( 'RootFxc1.sums', 'Test04\RootFxc1.sums' )
copyfile( 'RootFyc1.sums', 'Test04\RootFyc1.sums' )
copyfile( 'RootFzc1.sums', 'Test04\RootFzc1.sums' )
FileComp( 'Test04\DLC2.3_1.rflo', 'TestFiles\Test04\DLC2.3_1.rflo', 5, 'CertTest.comp' )
FileComp( 'Test04\RootFxc1.sums', 'TestFiles\Test04\RootFxc1.sums', 5, 'CertTest.comp' )
FileComp( 'Test04\RootFyc1.sums', 'TestFiles\Test04\RootFyc1.sums', 5, 'CertTest.comp' )
FileComp( 'Test04\RootFzc1.sums', 'TestFiles\Test04\RootFzc1.sums', 5, 'CertTest.comp' )

fprintf( '\n=======\nTest_05\n=======\n' );
MCrunch 'Test_05.mcru'
if ( ~exist( 'Test05', 'dir' ) )
   mkdir( 'Test05' )
end % if
copyfile( 'DLC2.3_1.extr'  , 'Test05\DLC2.3_1.extr'   )
copyfile( 'DLC2.3_1.pdfs'  , 'Test05\DLC2.3_1.pdfs'   )
copyfile( 'DLC2.3_1.rflo'  , 'Test05\DLC2.3_1.rflo'   )
copyfile( 'DLC2.3_1.stat'  , 'Test05\DLC2.3_1.stat'   )
copyfile( 'DLC2.3_2.extr'  , 'Test05\DLC2.3_2.extr'   )
copyfile( 'DLC2.3_2.pdfs'  , 'Test05\DLC2.3_2.pdfs'   )
copyfile( 'DLC2.3_2.rflo'  , 'Test05\DLC2.3_2.rflo'   )
copyfile( 'DLC2.3_2.stat'  , 'Test05\DLC2.3_2.stat'   )
copyfile( 'Test05_Agg.extr', 'Test05\Test05_Agg.extr' )
copyfile( 'Test05_Agg.pdfs', 'Test05\Test05_Agg.pdfs' )
copyfile( 'Test05_Agg.rflo', 'Test05\Test05_Agg.rflo' )
copyfile( 'Test05_Agg.stat', 'Test05\Test05_Agg.stat' )
FileComp( 'Test05\DLC2.3_1.extr'  , 'TestFiles\Test05\DLC2.3_1.extr'  , 5, 'CertTest.comp' )
FileComp( 'Test05\DLC2.3_1.pdfs'  , 'TestFiles\Test05\DLC2.3_1.pdfs'  , 5, 'CertTest.comp' )
FileComp( 'Test05\DLC2.3_1.rflo'  , 'TestFiles\Test05\DLC2.3_1.rflo'  , 5, 'CertTest.comp' )
FileComp( 'Test05\DLC2.3_1.stat'  , 'TestFiles\Test05\DLC2.3_1.stat'  , 5, 'CertTest.comp' )
FileComp( 'Test05\DLC2.3_2.extr'  , 'TestFiles\Test05\DLC2.3_2.extr'  , 5, 'CertTest.comp' )
FileComp( 'Test05\DLC2.3_2.pdfs'  , 'TestFiles\Test05\DLC2.3_2.pdfs'  , 5, 'CertTest.comp' )
FileComp( 'Test05\DLC2.3_2.rflo'  , 'TestFiles\Test05\DLC2.3_2.rflo'  , 5, 'CertTest.comp' )
FileComp( 'Test05\DLC2.3_2.stat'  , 'TestFiles\Test05\DLC2.3_2.stat'  , 5, 'CertTest.comp' )
FileComp( 'Test05\Test05_Agg.extr', 'TestFiles\Test05\Test05_Agg.extr', 5, 'CertTest.comp' )
FileComp( 'Test05\Test05_Agg.pdfs', 'TestFiles\Test05\Test05_Agg.pdfs', 5, 'CertTest.comp' )
FileComp( 'Test05\Test05_Agg.rflo', 'TestFiles\Test05\Test05_Agg.rflo', 5, 'CertTest.comp' )
FileComp( 'Test05\Test05_Agg.stat', 'TestFiles\Test05\Test05_Agg.stat', 5, 'CertTest.comp' )

fprintf( '\n=======\nTest_06\n=======\n' );
MCrunch 'Test_06.mcru'
if ( ~exist( 'Test06', 'dir' ) )
   mkdir( 'Test06' )
end % if
copyfile( 'Test_06.bins', 'Test06\Test_06.bins' )
copyfile( 'Test_06.psds', 'Test06\Test_06.psds' )
FileComp( 'Test06\Test_06.bins', 'TestFiles\Test06\Test_06.bins', 5, 'CertTest.comp' )
FileComp( 'Test06\Test_06.psds', 'TestFiles\Test06\Test_06.psds', 5, 'CertTest.comp' )

fprintf( '\n=======\nTest_07\n=======\n' );
MCrunch 'Test_07.mcru'
if ( ~exist( 'Test07', 'dir' ) )
   mkdir( 'Test07' )
end % if
copyfile( 'Test_07_Agg.stat'    , 'Test07\Test_07_Agg.stat'     )
copyfile( 'Test_07_Agg.life'    , 'Test07\Test_07_Agg.life'     )
copyfile( 'DLC1.1_01_small.stat', 'Test07\DLC1.1_01_small.stat' )
copyfile( 'DLC1.1_07_small.stat', 'Test07\DLC1.1_07_small.stat' )
copyfile( 'DLC1.1_13_small.stat', 'Test07\DLC1.1_13_small.stat' )
copyfile( 'DLC1.1_19_small.stat', 'Test07\DLC1.1_19_small.stat' )
copyfile( 'DLC1.1_25_small.stat', 'Test07\DLC1.1_25_small.stat' )
copyfile( 'DLC1.1_31_small.stat', 'Test07\DLC1.1_31_small.stat' )
copyfile( 'DLC1.1_37_small.stat', 'Test07\DLC1.1_37_small.stat' )
copyfile( 'DLC1.1_43_small.stat', 'Test07\DLC1.1_43_small.stat' )
copyfile( 'DLC1.1_49_small.stat', 'Test07\DLC1.1_49_small.stat' )
FileComp( 'Test07\Test_07_Agg.stat'    , 'TestFiles\Test07\Test_07_Agg.stat'    , 5, 'CertTest.comp' )
FileComp( 'Test07\Test_07_Agg.life'    , 'TestFiles\Test07\Test_07_Agg.life'    , 5, 'CertTest.comp' )
FileComp( 'Test07\DLC1.1_01_small.stat', 'TestFiles\Test07\DLC1.1_01_small.stat', 5, 'CertTest.comp' )
FileComp( 'Test07\DLC1.1_07_small.stat', 'TestFiles\Test07\DLC1.1_07_small.stat', 5, 'CertTest.comp' )
FileComp( 'Test07\DLC1.1_13_small.stat', 'TestFiles\Test07\DLC1.1_13_small.stat', 5, 'CertTest.comp' )
FileComp( 'Test07\DLC1.1_19_small.stat', 'TestFiles\Test07\DLC1.1_19_small.stat', 5, 'CertTest.comp' )
FileComp( 'Test07\DLC1.1_25_small.stat', 'TestFiles\Test07\DLC1.1_25_small.stat', 5, 'CertTest.comp' )
FileComp( 'Test07\DLC1.1_31_small.stat', 'TestFiles\Test07\DLC1.1_31_small.stat', 5, 'CertTest.comp' )
FileComp( 'Test07\DLC1.1_37_small.stat', 'TestFiles\Test07\DLC1.1_37_small.stat', 5, 'CertTest.comp' )
FileComp( 'Test07\DLC1.1_43_small.stat', 'TestFiles\Test07\DLC1.1_43_small.stat', 5, 'CertTest.comp' )
FileComp( 'Test07\DLC1.1_49_small.stat', 'TestFiles\Test07\DLC1.1_49_small.stat', 5, 'CertTest.comp' )

edit( 'CertTest.comp' )

beep;
