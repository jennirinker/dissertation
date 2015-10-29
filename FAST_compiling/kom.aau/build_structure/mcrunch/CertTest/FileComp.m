function FileComp( File1, File2, MaxDiff, DiffFile )
% Compare File1 to File2 and add differling lines to DiffFile.
% Its primary purpose is to document differences in two versions of MCrunch.
%
% Syntax is:  FileComp( File1, File2, MaxDiff, DiffFile );
%
% where:
%     File1    - The first file to compare.
%     File2    - The second file to compare.
%     MaxDiff  - The maximum number of differences to document in DiffFile.
%     DiffFile - The file that has the differences appended to it.
%
% Example:
%     FileComp( 'DLC2.3_1.extr', 'TestFiles\DLC2.3_1.extr', 5, 'CertTest.comp' );
%
% See also CertTest, MCrunch


      % Open files.

   Un1 = fopen( File1, 'rt' );
   if ( Un1 < 0 )
      beep
      error( '  Could not open "%s" for reading.', File1 );
   end

   Un2 = fopen( File2, 'rt' );
   if ( Un2 < 0 )
      beep
      error( '  Could not open "%s" for reading.', File2 );
   end

   UnD = fopen( DiffFile, 'at' );
   if ( UnD < 0 )
      beep
      error( '  Could not open "%s" for appending.', DiffFile );
   end


      % Compare the files.

   fprintf( UnD, '\n==================================================================================\n' );
   fprintf( UnD, 'File differences between "%s" and "%s":\n', File1, File2 );

   Line1   = fgets( Un1 );
   Line2   = fgets( Un2 );
   Line    = 1;
   NumDiff = 0;

   while ( ischar( Line1 ) && ischar( Line2 ) )

      if ( ~strcmp( Line1, Line2 ) )

         fprintf( UnD, '\n  Line #%d\n  1: %s  2: %s', Line, Line1, Line2 );
         NumDiff = NumDiff + 1;

         if ( NumDiff > MaxDiff )
            fprintf( UnD, '\n  Maximum number of differences reached.\n' );
            break
         end % if

      end % if ( ~strcmp( Line1, Line2 ) )

      Line1 = fgets( Un1 );
      Line2 = fgets( Un2 );
      Line  = Line + 1;

   end % while


      % Close the files.

   fclose( Un1 );
   fclose( Un2 );
   fclose( UnD );


   return

end % function FileComp( File1, File2, MaxDiff, DiffFile )