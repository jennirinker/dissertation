#!/bin/perl -w

# createpage.pl - Creates the TurbSim web page.

use strict;


# The following variables should be updated for every new version of the code.

my $arch_file     = "TurbSim_v1.06.00.exe";
my $arch_ver      = "v1.06.00";
my $compiler      = "Intel&reg; Visual Fortran Compiler XE 12.1.3.300 [IA-32]";
my $compiler64    = "Intel&reg; Visual Fortran Compiler XE 12.1.3.300 [Intel&reg; 64]";

my $guide_ver     = "v1.06.00";                            # version number corresponding to the TurbSim user's guide
my $over_ver      = "v1.21";                               # version number corresponding to the TurbSim overview document
my $nwtclib_ver   = "v1.04.01";                            # the version number of the nwtc library


# The following variables should rarely change.

my $bin_file      = "TSM_structures.exe";
my $change_file   = "ChangeLog.txt";
my $dns_file      = "TSM_DNS_structures.exe";
my $email1        = '';
my $email2        = 'bonnie%2Ejonkman%40nrel%2Egov';
my $guide_file    = "TurbSim.pdf";
my $name1         = "Neil Kelley";
my $name2         = "Bonnie Jonkman";
my $nwtclib_path  = "http://wind.nrel.gov/designcodes/miscellaneous/nwtc_subs/";
my $over_file     = "TurbSimOverview.pdf";
my $path          = "Y:/Wind/WindWeb/designcodes/preprocessors/turbsim/";
my $program       = "TurbSim";
my $web_page      = "index.html";
my $web_path      = "http://wind.nrel.gov/designcodes/preprocessors/turbsim/";
my $pdfIcon       = "http://wind.nrel.gov/designcodes/icon_pdf.gif";


# Set file permissions and ownership.

chmod 0664, $arch_file, $bin_file, $dns_file;
chmod 0664, $guide_file, $change_file, $web_page, $over_file;
chmod 0774, $0;


# Generate file size information.
my $change_size = file_sz($path.$change_file);  #change log
my $over_size   = file_sz($path.$over_file);    #Overview
my $guide_size  = file_sz($path.$guide_file);   #user's guide
my $arch_size   = file_sz($path.$arch_file);    #archive
my $bin_size    = file_sz($path.$bin_file);     #cohverent structures


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

<title>NWTC Computer-Aided Engineering Tools ($program)</title>

<!--#include virtual="/header.html"-->

</head>

<body>

<!--#include virtual="/designcodes/title.html"-->

<table border="0" cellpadding="5" width="100%">

  <tr align="left" valign="top">
     <th>Contents&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</th>
     <th>&nbsp</th>
     <th>$program</th></tr>
  <tr align="left" valign="top">
     <td><!--#include virtual="/designcodes/content_links_all.html"--></td>
     <td>&nbsp;</td>
     <td class="small">

  <i>A stochastic, full-field, turbulent-wind simulator<BR>
   for use with InflowWind/AeroDyn-based simulation tools</i>

  <br><br>
  <strong> by $name1
  and <a href="mailto:$email2">$name2</a></strong>
  <br>National Wind Technology Center


  <p>The $program stochastic inflow turbulence tool has been developed to provide a numerical
  simulation of a full-field flow that contains coherent turbulence structures that reflect the
  proper spatiotemporal turbulent velocity field relationships seen in instabilities associated with
  nocturnal boundary layer flows and which are not represented well by the IEC Normal Turbulence
  Models (NTM).  Its purpose is to provide the wind turbine designer with the ability to drive
  design code (e.g., FAST or MSC.Adams&reg) simulations of advanced turbine designs with simulated inflow turbulence
  environments that incorporate many of the important fluid dynamic features known to adversely affect
  turbine aeroelastic response and loading.</p>

  <p>$program supports all of the features found in the previous SNLWIND-3D<sup>1</sup> and SNwind<sup>2</sup> inflow
  turbulence simulation tools that had their roots in Paul Veers&#39; original stochastic wind simulator, SNLWIND<sup>3</sup>.
  All of the spectral models available in SNLWIND-3D (SMOOTH, WF_UPW, WF_07D, and WF_14D) as well as those available
  in SNwind (the IEC Kaimal and von Karman Normal Turbulence Models (NTM)) are included in $program.  In addition,
  two new spectral models, NWTCUP and GP_LLJ, have been added to $program.  The NWTCUP model is based on conditions
  observed at the National Wind Technology Center (NWTC); the GP_LLJ model is based on conditions associated with Great Plains
  low-level jet streams.  $program is more efficient than its predecessors in terms of both CPU and memory usage.</p>


  <p>$program provides the ability to efficiently generate randomized coherent turbulent structures
  that are superimposed on the more random background turbulent field as produced by one of the diabatic
  (non-neutral) spectral models; i.e., GP_LLJ, NWTCUP, SMOOTH, WF_UPW, WF_07D, or WF_14D. Modeled as a combination of
  non-homogenous Poisson and Lognormal Stochastic Processes, the randomized scaling of the coherent
  structures in this version of $program are based on measurements derived from the Lamar 120-m tower used as part of the
  Lamar Low-Level Jet Project (LLLJP) when the GP_LLJ or SMOOTH models are specified.
  The NWTCUP model uses scaling derived from measurements taken from the 5-sonic anemometer
  array upwind of the NWTC ART Turbine as part of the Long-Term Inflow and Structural Testing (LIST) Program, and
  the three wind farm models (WF_UPW, WF_07D, and WF_14D) use scaling derived from measurements taken from upwind, within, and
  downwind of a large wind farm in San Gorgonio Pass, California.
  </p>

  <p> In version 1.06.00, $program was extended to model turbulence in water. The TIDAL spectral model, which simulates
  turbulence in a tidal channel, is based on measurements taken near Marrowstone Island in Puget Sound Washington. </p>

  <ol>
  <li> Kelley, N.D. (November 1992). <i>Full Vector (3-D) Inflow Simulation in Natural and Wind Farm Environments
           Using an Expanded Version of the SNLWIND (Veers) Turbulence Code</i>. NREL/TP-442-5225.
           Golden, CO: National Renewable Energy Laboratory.
  <li> Buhl, M.L. Jr. (October 2001). <i>SNwind User&#39;s Guide</i>. NREL/EL-500-30121.
           Golden, CO: National Renewable Energy Laboratory.
  <li> Veers, P.S. (March 1988). <i>Three-Dimensional Wind Simulation</i>, SAND88-0152.
           Albuquerque, NM: Sandia National Laboratories.
  </ol>



  <p><b>You may download the following files from our server:</b></p>
  <!--#config timefmt="%e-%B-%Y" -->

  <ul>

    <li><a href="$change_file">$program Change Log</a> ($arch_ver, $change_size, <!--#flastmod virtual="$change_file" -->)

    <p>This is a list of changes made to the code. Look at this text
    file to see if we&#39;ve made worthwhile changes since you received your
    previous version of $program.</p>


    <li><a href="$over_file">$program Overview<img src="$pdfIcon" style="border:0"
    width="13" height="14" class="arrowicon" alt="PDF"></a> ($over_ver, $over_size, <!--#flastmod virtual="$guide_file" -->)

    <p>This is an overview of $program features. Please refer to it when trying to
    understand what the program does.</p>


    <li><a href="$guide_file">$program User&#39;s Guide<img src="$pdfIcon" style="border:0" width="13"
    height="14" class="arrowicon" alt="PDF"></a> ($guide_ver, $guide_size, <!--#flastmod virtual="$guide_file" -->)

    <p>This is the $program User&#39;s Guide. Please refer to it when trying to
    understand how the program works.</p>


    <li><a href="${web_path}downloaders/$arch_file">$program Archive (EXE $arch_size)</a> ($arch_ver, $arch_size, <!--#flastmod virtual="$arch_file" -->)

    <p>This is a self-extracting archive of $program. It contains the
    $program executable file, sample input file, installation-verification
    test procedure, change log, source code, and user&#39;s guide. It runs on
    all 32-bit Windows&reg; platforms. We created the executable file with
    $compiler.  The archive also includes a $program executable file for 64-bit
    Windows&reg; platforms.  It was created with $compiler64. The code is compiled with
    <a href="$nwtclib_path">NWTC Subroutine Library</a> $nwtclib_ver.</p>

    <li><a href="${web_path}downloaders/$bin_file">$program Coherent Structures Archive (EXE $bin_size)</a> (<!--#flastmod virtual="$bin_file" -->)

    <p>This is a self-extracting archive of $program coherent structures. It contains
    the event definition files for both LES and DNS structure types
    as well as the associated binary data that is used to add coherent structures to
    background turbulence in AeroDyn. It runs on all 32-bit Windows&reg; platforms.
    You will need these files in order to generate coherent structures in $program.</p>


  </ul>

  <p><b>If you want to refer to the $program website in a report, here is a reference
  you can use:</b></p>

  <!--#config timefmt="%e-%B-%Y" -->
  <blockquote><p><i>NWTC Computer-Aided Engineering Tools ($program by $name1, $name2)</i>.
  $web_path. Last modified <!--#flastmod virtual="$web_page" -->; accessed <!--#echo var="DATE_LOCAL" -->. </p></blockquote>

  <!--#config timefmt="%A, %e-%B-%Y at %T %Z" -->
  <br><br><hr><i class="s70">This page was generated by CreatePage.pl on $date_time. <br> Last modification: <!--#flastmod virtual="$web_page" --> </i>

  <!--#include virtual="/validator.html"-->

  </td></tr>
</table>

<p></p><!--#include virtual="/nrel_footer.html"-->

</body>
</html>
EOF
# ============================= End Here document =======================================

close INDEX;

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

#*******************************************************************
#*******************************************************************
# This routine returns the file size.

sub file_sz
{

   my( $fileName ) = @_;
   my $fs = "0 B";


   if ( -e $fileName ) {
      my @dummy  = stat $fileName;
      $fs     = $dummy[7]/1000;

      if ( $fs > 1000 ) {
         $fs = sprintf( "%.1f MB", $fs/1000 );
      } else {
         $fs = sprintf( "%.0f KB", $fs );
      }

      # EERE says we're supposed to remove trailing zeros, too...

   } else {
      print ("\n  WARNING: $fileName does not exist.\n");
   }

   return $fs;
}

__END__

