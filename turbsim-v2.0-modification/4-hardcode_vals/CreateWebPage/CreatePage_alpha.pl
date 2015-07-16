#!/bin/perl -w

# createpage.pl - Creates the TurbSim web page.

use strict;


# The following variables should be updated for every new version of the code.

my $arch_file   = "TurbSim_v1.06.00a-bjj.exe";
my $arch_ver    = "v1.06.00a-bjj";
my $compiler    = "Intel&reg; Visual Fortran v10.1.024"; #9.0.030";
my $compiler64  = "Intel&reg; Visual Fortran 64 Compiler XE 12.1.0.233"; 
my $guide_ver   = "v1.50";
my $over_ver    = "v1.21";


# The following variables should rarely change.

my $change_file = "ChangeLog.txt";
my $email1      = 'neil%2Ekelley%40nrel%2Egov';
my $email2      = 'bonnie%2Ejonkman%40nrel%2Egov';
my $guide_file  = "TurbSim.pdf";
my $name1       = "Neil Kelley";
my $name2       = "Bonnie Jonkman";
my $over_file   = "TurbSimOverview.pdf";
my $path        = "Y:/Wind/WindWeb/designcodes/preprocessors/turbsim/alpha/";
my $program     = "TurbSim";
my $web_page    = "index.html";
my $web_path    = "http://wind.nrel.gov/designcodes/preprocessors/turbsim/alpha/";
my $web_path_b  = "http://wind.nrel.gov/designcodes/preprocessors/turbsim/";
my $pdfIcon     = "http://wind.nrel.gov/designcodes/icon_pdf.gif";


# Set file permissions and ownership.

#chmod 0664, $arch_file, $bin_file, $over_file;
chmod 0664, $arch_file;
chmod 0664, $guide_file, $change_file, $web_page;
chmod 0774, $0;


# Generate change-log file information.

my @dummy       = stat $path.$change_file;
my $change_size = $dummy[7]/1024;

if ( $change_size > 1024 ) {
	$change_size = sprintf( "%.1f MB", $change_size/1024 );
} else {
	$change_size = sprintf( "%.0f KB", $change_size );
}

my $change_date = &date_str( gmtime( $dummy[9] ), "long" );


# Generate User's Guide file information.

   @dummy      = stat $path.$guide_file;
my $guide_size = $dummy[7]/1024;

if ( $guide_size > 1024 ) {
	$guide_size = sprintf( "%.1f MB", $guide_size/1024 );
} else {
	$guide_size = sprintf( "%.0f KB", $guide_size );
}

my $guide_date = &date_str( gmtime( $dummy[9] ), "long" );


# Generate archive file information.

   @dummy     = stat $path.$arch_file;
my $arch_size = $dummy[7]/1024;

if ( $arch_size > 1024 ) {
	$arch_size = sprintf( "%.1f MB", $arch_size/1024 );
} else {
	$arch_size = sprintf( "%.0f KB", $arch_size );
}

my $arch_date = &date_str( gmtime( $dummy[9] ), "long" );


# Generate date and time.

my $date      = date_now();
my $date_time = sprintf "%s at %s", $date, time_now();


# Create the web page.

open INDEX, ">$path$web_page" or die "Can't open $path$web_page for writing.";

# ============================ Begin Here document ======================================
print  INDEX <<EOF;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">

<html>

<head>

<title>NWTC Computer-Aided Engineering Tools ($program Alpha Testing)</title>

<!--#include virtual="/header.html"-->

</head>

<body>

<!--#include virtual="/designcodes/title.html"-->

<table border="0" cellpadding="5" width="100%">

  <tr align="left" valign="top"><th>Contents&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</th>
                                <th>&nbsp</th><th>$program Alpha-Testers Page</th></tr>
  <tr align="left" valign="top">
  <td><!--#include virtual="/link_home.html"-->
  <br><!--#include virtual="/link_contents_01.html"-->
  <table border="0" cellpadding="5" width="90%">
    <tr align="left" valign="top">
    <td><!--#include virtual="/designcodes/link_all_CAE_tools.html"-->
  </table>
      <!--#include virtual="/link_contents_02.html"-->
      <!--#include virtual="/link_contents_03.html"-->
      <!--#include virtual="/link_contents_04.html"-->
      <!--#include virtual="/link_contents_05.html"-->
      <!--#include virtual="/link_contents_06.html"-->
      <!--#include virtual="/link_contents_07.html"-->
      <!--#include virtual="/link_contents_08.html"-->
      <!--#include virtual="/link_contents_09.html"-->
      <!--#include virtual="/link_contents_10.html"-->
  <br><!--#include virtual="/link_other.html"-->
  </td>
  <td>
  &nbsp;
  </td>
  <td class="small">  

  <i>A stochastic, full-field, turbulent-wind simulator<BR>
   for use with InflowWind/AeroDyn-based CAE tools</i>

  <br><br>
  <strong> by $name1
  and <a href="mailto:$email2">$name2</a></strong>
  <br>National Wind Technology Center

  <p> This page is for alpha testers of $program.  It contains links to an alpha version of the code, which <i>has
      not been fully tested</i>.  The latest beta version of the code and additional code documentation 
      can be found on the main <a href="$web_path_b">
      $program webpage</a>. </p>

  <p> Please send comments regarding your experience with this alpha version to <a href="mailto:$email2">$name2</a>. </p>

  <p><b>Alpha testers may download the following files from our server:</b></p>

  <ul>

    <li><a href="$change_file">$program Change Log</a> ($arch_ver, $change_size, $change_date)

    <p>This is a list of changes made to the code.&nbsp; Look at this text
    file to see what changes have been made since you received your
    previous version of $program.</p>

    <li><a href="downloaders/$arch_file">$program Archive (EXE $arch_size)</a> ($arch_ver, $arch_size, $arch_date)

    <p>This is a self-extracting archive of an alpha version of $program.&nbsp; It contains the
    $program executable file, sample input file, installation-verification
    test procedure, change log, source code, and user's guide.&nbsp; It runs on
    all 32-bit Windows&reg; platforms.&nbsp; We created the executable file with
    $compiler. The archive also includes a $program executable file for 64-bit
    Windows&reg; platforms.  It was created with $compiler64.</p>

  </ul>


  <br><br><hr><i class="s70">This page was generated by createpage.pl on $date_time.</i>

  <p>
    <a href="http://validator.w3.org/check?uri=$web_path$web_page">
      <img style="border:0;width:88px;height:31px"
         src="http://www.w3.org/Icons/valid-html401"
         alt="Valid HTML 4.01!" height="31" width="88">
    </a>
    <br>
    <a href="http://jigsaw.w3.org/css-validator/validator?uri=$web_path$web_page">
      <img style="border:0;width:88px;height:31px"
         src="http://jigsaw.w3.org/css-validator/images/vcss"
         alt="Valid CSS!">
    </a>
  </p>

  </td></tr>
</table>

<p></p><!--#include virtual="/nrel_footer.html"-->

</body>
</html>
EOF
# ============================= End Here document =======================================

close INDEX;

#    <li><a href="$guide_file">$program User's Guide</a> ($guide_ver, $guide_size, $guide_date)
#
#    <p>This is the $program User's Guide.&nbsp; Please refer to it when trying to
#    understand how the program works.</p>
#

# End of main function.

#*******************************************************************
#*******************************************************************
# This routine returns the current date in the form "dd-Mon-ccyy" or "dd-Month-ccyy".

sub date_now
{
   return date_str( localtime, "long" );
}


#*******************************************************************
#*******************************************************************
# This routine converts a date ($mday,$mon,$year) format and converts it to the form "dd-Mon-ccyy" or "dd-Month-ccyy".

sub date_str
{
   my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst, $long) = @_;

   my $mon_str;

   if ( defined $long and $long eq "long" ) {
      $mon_str = ("January","February","March","April","May","June","July","August","September","October","November","December")[$mon];
	} else {
      $mon_str = ("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")[$mon];
	}

   return sprintf( "%2.2d-%s-%4.4d", $mday, $mon_str, $year+1900 );
}


#*******************************************************************
#*******************************************************************
# This routine returns the time in the form "hh:mm:ss".

sub time_now
{
   my( $sec, $min, $hour, $mday, $mon, $year ) = localtime;

   return sprintf( "%2.2d:%2.2d:%2.2d", $hour, $min, $sec );
}

__END__

