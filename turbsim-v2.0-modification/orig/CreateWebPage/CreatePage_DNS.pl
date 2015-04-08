#!/bin/perl -w

# createpage.pl - Creates the TurbSim web page.

use strict;


# The following variables should be updated for every new version of the code.

my $arch_file   = "TurbSim_v110.exe";
my $arch_ver    = "v1.10";
my $compiler    = "Intel&reg; Visual Fortran v9.0.030";
my $guide_ver   = "v1.10";
my $over_ver    = "v1.10";


# The following variables should rarely change.

my $bin_file    = "TSM_structures.exe"; 
my $change_file = "ChangeLog.txt";
my $dns_file    = "TSM_DNS_structures.exe";
my $email1      = 'neil%2Ekelley%40nrel%2Egov';
my $email2      = 'bonnie%2Ejonkman%40nrel%2Egov';
my $guide_file  = "TurbSim.pdf";
my $name1       = "Neil Kelley";
my $name2       = "Bonnie Jonkman";
my $over_file   = "TurbSimOverview.pdf";
my $path        = "Y:/Wind/WindWeb/designcodes/preprocessors/turbsim/";
my $program     = "TurbSim";
my $web_page    = "index.html";
my $web_path    = "http://wind.nrel.gov/designcodes/preprocessors/turbsim/";
my $pdfIcon     = "http://wind.nrel.gov/designcodes/icon_pdf.gif";

# Set file permissions and ownership.

chmod 0664, $arch_file, $bin_file, $dns_file;
chmod 0664, $guide_file, $change_file, $web_page, $over_file;
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


# Generate Overview file information.

   @dummy      = stat $path.$over_file;
my $over_size = $dummy[7]/1024;

if ( $over_size > 1024 ) {
	$over_size = sprintf( "%.1f MB", $over_size/1024 );
} else {
	$over_size = sprintf( "%.0f KB", $over_size );
}

my $over_date = &date_str( gmtime( $dummy[9] ), "long" );


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


# Generate coherent structure archive file information.

   @dummy     = stat $path.$bin_file;
my $bin_size = $dummy[7]/1024;

if ( $bin_size > 1024 ) {
	$bin_size = sprintf( "%.1f MB", $bin_size/1024 );
} else {
	$bin_size = sprintf( "%.0f KB", $bin_size );
}

my $bin_date = &date_str( gmtime( $dummy[9] ), "long" );


# Generate DNS coherent structure archive file information.

   @dummy     = stat $path.$dns_file;
my $dns_size = $dummy[7]/1024;

if ( $dns_size > 1024 ) {
	$dns_size = sprintf( "%.1f MB", $dns_size/1024 );
} else {
	$dns_size = sprintf( "%.0f KB", $dns_size );
}

my $dns_date = &date_str( gmtime( $dummy[9] ), "long" ); 


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

<title>NWTC Design Codes ($program)</title>

<!--#include virtual="/header.html"-->

</head>

<body>

<!--#include virtual="/designcodes/title.html"-->

<table border="0" cellpadding="5" width="100%">

  <tr align="left" valign="top"><th>Contents&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</th>
                                <th>&nbsp</th><th>$program</th></tr>
  <tr align="left" valign="top">
  <td><!--#include virtual="/link_home.html"-->
  <br><!--#include virtual="/link_amestest.html"-->
  <br><!--#include virtual="/link_annex-xx.html"-->
  <br><!--#include virtual="/link_cert_stds.html"-->
  <br><!--#include virtual="/link_designcodes.html"-->
  <table border="0" cellpadding="5" width="90%">
    <tr align="left" valign="top">
    <td><!--#include virtual="/designcodes/link_policy.html"-->
    <br><!--#include virtual="/designcodes/link_preprocessors.html"-->
    <br><!--#include virtual="/designcodes/link_simulators.html"-->
    <br><!--#include virtual="/designcodes/link_postprocessors.html"-->
    <br><!--#include virtual="/designcodes/link_miscellaneous.html"-->
    <br><!--#include virtual="/designcodes/link_setup.html"-->
    <br><!--#include virtual="/designcodes/link_advice.html"-->
    <br><!--#include virtual="/designcodes/link_papers.html"-->
  </table>
      <!--#include virtual="/link_dyno.html"-->
  <br><!--#include virtual="/link_furling.html"-->
  <br><!--#include virtual="/link_mt15.html"-->
  <br><!--#include virtual="/link_osu.html"-->
  <br><!--#include virtual="/link_swrt.html"-->
  <br><!--#include virtual="/link_other.html"-->
  </td>
  <td>
  &nbsp;
  </td>
  <td class="small">

  <i>A stochastic, full-field, turbulent-wind simulator for use<BR>
  with the AeroDyn-based design
  codes <BR>(YawDyn, FAST, and MSC.ADAMS&reg;)</i>

  <br><br>
  <strong> by <a href="mailto:$email1">$name1</a>
  and <a href="mailto:$email2">$name2</a></strong>
  <br>National Wind Technology Center


  <p>The $program stochastic inflow turbulence code has been developed to provide a numerical 
  simulation of a full-field flow that contains coherent turbulence structures that reflect the 
  proper spatiotemporal turbulent velocity field relationships seen in instabilities associated with 
  nocturnal boundary layer flows and which are not represented well by the IEC Normal Turbulence 
  Models (NTM).  Its purpose is to provide the wind turbine designer with the ability to drive 
  design code (YawDyn, FAST, or MSC.ADAMS&reg) simulations of advanced turbine designs with simulated inflow turbulence
  environments that incorporate many of the important fluid dynamic features known to adversely affect
  turbine aeroelastic response and loading.</p>

  <p>The $program code supports all of the features found in the previous SNLWIND-3D<sup>1</sup> and SNwind<sup>2</sup> inflow
  turbulence simulator codes that had their roots in Paul Veers' original stochastic wind simulator SNLWIND<sup>3</sup>.
  All of the spectral models available in SNLWIND-3D (SMOOTH, WF_UPW, WF_07D, and WF_14D) as well as the 
  IEC Kaimal and von Karman Normal Turbulence Models (NTM) that are provided in SNwind are included.  
  Two new spectral models, NWTCUP and GP_LLJ, have been added to the TurbSim code.  These models are based on conditions 
  observed the National Wind Technology Center (NWTC) and conditions associated with Great Plains low-level jet streams
  respectively.  $program is more efficient than its predecessors in terms of both CPU and memory usage.</p>
  

  <p>The $program code provides the ability to efficiently generate randomized coherent turbulent structures 
  that are superimposed on the more random background turbulent field as produced by one of the diabatic 
  (non-neutral) spectral models; i.e., GP_LLJ, NWTCUP, SMOOTH, WF_UPW, WF_07D, or WF_14D. Modeled as a combination of 
  non-homogenous Poisson and Lognormal Stochastic Processes, the randomized scaling of the coherent 
  structures in this version of TurbSim are based on measurements derived from the Lamar 120-m tower used as part of the
  Lamar Low-Level Jet Project (LLLJP) when the GP_LLJ model is specified. 
  The other models use scaling derived from measurements taken from the 5-sonic anemometer 
  array upwind of the NWTC ART Turbine as part of the Long-Term Inflow and Structural Testing (LIST) Program. 
  </p>


  <ol>  
  <li> Kelley, N.D. (November 1992). <i>Full Vector (3-D) Inflow Simulation in Natural and Wind Farm Environments
           Using an Expanded Version of the SNLWIND (Veers) Turbulence Code</i>. NREL/TP-442-5225. 
           Golden, CO: National Renewable Energy Laboratory. 
  <li> Buhl, M.L. Jr. (October 2001). <i>SNwind User's Guide</i>. NREL/EL-500-30121. 
           Golden, CO: National Renewable Energy Laboratory.
  <li> Veers, P.S. (March 1988). <i>Three-Dimensional Wind Simulation</i>, SAND88-0152. 
           Albuquerque, NM: Sandia National Laboratories. 
  </ol>



  <p><b>You may download the following files from our server:</b></p>

  <ul>

    <li><a href="$change_file">$program Change Log</a> ($arch_ver, $change_size, $change_date)

    <p>This is a list of changes made to the code.&nbsp; Look at this text
    file to see if we've made worthwhile changes since you received your
    previous version of $program.</p>

    <li><a href="$over_file">$program Overview<img src="$pdfIcon" style="border:0" width="13" 
    height="14" class="arrowicon" alt="PDF"></a></a> ($over_ver, $over_size, $over_date)

    <p>This is an overview of $program features.&nbsp; Please refer to it when trying to
    understand what the program does and how to use it.</p>

    <li><a href="$guide_file">$program User's Guide<img src="$pdfIcon" style="border:0" 
    width="13" height="14" class="arrowicon" alt="PDF"></a></a> ($guide_ver, $guide_size, $guide_date)

    <p>This is the $program User's Guide.&nbsp; Please refer to it when trying to
    understand how the program works.</p>


    <li><a href="$arch_file">$program Archive (EXE $arch_size)</a> ($arch_ver, $arch_size, $arch_date)

    <p>This is a self-extracting archive of $program.&nbsp; It  contains the
    $program executable file, sample input file, installation-verification
    test procedure, change log, source code, and user's guide.&nbsp; It runs on
    all 32-bit Windows&reg; platforms.&nbsp; We created the executable file with
    $compiler.</p>

    <li><a href="$bin_file">$program Coherent Structures Archive (EXE $bin_size)</a> ($bin_date)

    <p>This is a self-extracting archive of $program coherent structures.&nbsp; It contains 
    the event definition files for <i>both</i> LES and DNS structure types
    as well as the associated binary data that is used to add coherent structures to 
    background turbulence in AeroDyn.&nbsp; It runs on all 32-bit Windows&reg; platforms.
    You will need these files in order to generate coherent structures in $program.</p>

    <li><a href="$dns_file">$program DNS Coherent Structures Archive (EXE $dns_size)</a> ($dns_date)

    <p>This is a self-extracting archive of $program DNS coherent structures.&nbsp; This archive
    is a subset of the $program Coherent Structures Archive, and is intended for people who have 
    already downloaded and installed the LES files that were distributed with $program v1.00a.&nbsp; 
    The archive contains event definition files and associated binary data for <i>only</i> 
    the DNS structure types.&nbsp; </p>

  </ul>

  <p><b>If you want to refer to the $program website in a report, here is a reference
  you can use:</b></p>

  <blockquote><p><i>NWTC Design Codes ($program by $name1, $name2)</i>.&nbsp;
  $web_path.&nbsp; Last modified $date; accessed $date. </p></blockquote>


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

