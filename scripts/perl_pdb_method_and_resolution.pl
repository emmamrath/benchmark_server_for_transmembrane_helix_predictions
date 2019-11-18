#!/usr/bin/perl   -w

# perl perl_pdb_method_and_resolution.pl -infile opm_feb2012.unique_real_seq_opm_in_pdb.pdbid

# Read in a list of PDBIDs. 
# For each PDBID, read info from pdb : http://www.rcsb.org/pdb/explore/explore.do?structureId=1zoy
# Output the method (x-ray diffraction or nmr) and resolution (x-ray diffraction only).

# input file :
# 1i78_A
# 1m56_B

# output file :
# 1zoy_A:::::METHOD,,,,,X-RAY DIFFRACTION;;;;;RESOLUTION,,,,,2.40;;;;;:::::
# 1zoy_B:::::METHOD,,,,,X-RAY DIFFRACTION;;;;;RESOLUTION,,,,,2.40;;;;;:::::
# 2lbv_A:::::METHOD,,,,,SOLUTION NMR;;;;;:::::

##################################################

use strict;
use Getopt::Long;
use LWP;
use HTTP::Request::Common;
use Data::Dumper;



##################################################

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
$year += 1900;
$mon += 1;
$mon = sprintf("%02d", $mon);
$mday = sprintf("%02d", $mday);
$hour = sprintf("%02d", $hour);
$min = sprintf("%02d", $min);
$sec = sprintf("%02d", $sec);
my $time_str = "$year-$mon-$mday.$hour:$min:$sec";

# constants

my $create_outfile2 = 1;
# http://www.rcsb.org/pdb/explore/explore.do?structureId=1zoy
my $url_1_for_submit = 'http://www.rcsb.org/pdb/explore/explore.do?structureId=';

# global variables

my $max_num_retries_of_query1 = 10;
my $count_input_pdbids = 0;
my $output_info;
my $referer = 'perl program';

# get the input parameters

my $input_file;
GetOptions( "infile=s" => \$input_file );
if (!defined $input_file) {
	die "Usage $0 -infile INPUTFILE\n";
}

# Open output files

my $output_file = "$input_file.pdb_method_resolution";
my $output_file_2 = "$input_file.pdb_method_html";
my $redo_file = "$input_file.pdb_method_redo";

open( OUTFILE, ">$output_file") or
	die "Cannot open $output_file for writing\n";
if ($create_outfile2 == 1) {
	open( OUTFILE2, ">$output_file_2") or
		die "Cannot open $output_file_2 for writing\n";
}
open( REDOFILE, ">$redo_file") or
	die "Cannot open $redo_file for writing\n";

# read input file

open INFILE, $input_file or die $!;
my @input_lines = <INFILE>;

# process each input line

foreach my $input_line (@input_lines) {
	chomp $input_line;
	$input_line = trim($input_line);
	if ($input_line ne '') {
		my @bits = split(/:::::/, $input_line);
		my $pdbid_and_chain = $bits[0];
		my $pdbid = substr( $pdbid_and_chain, 0, 4 );
		process_1_pdbid( $pdbid, $pdbid_and_chain );
	}
}

# close the output files

#print OUTFILE "\n";
close OUTFILE;
if ($create_outfile2 == 1) {
	print OUTFILE2 "\n";
	close OUTFILE2;
}
close REDOFILE;


##################################################
sub process_1_pdbid {

	my $pdbid = shift;
	my $pdbid_and_chain = shift;

	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	$year += 1900;
	$mon += 1;
	$mon = sprintf("%02d", $mon);
	$mday = sprintf("%02d", $mday);
	$hour = sprintf("%02d", $hour);
	$min = sprintf("%02d", $min);
	$sec = sprintf("%02d", $sec);
	$time_str = "$year-$mon-$mday.$hour:$min:$sec";

	$count_input_pdbids++;
	print         "$time_str processing input number $count_input_pdbids, sequence $pdbid_and_chain...\n";
	if ($create_outfile2 == 1) {
		print OUTFILE2 "$time_str Processing $pdbid_and_chain\n\n";
	}

	my $ua = LWP::UserAgent->new;
	my $newagent = 'test search program (' . $ua->agent . ')';
	$ua->agent($newagent);

	my $url_for_submit = $url_1_for_submit . $pdbid;

	my $request = GET ( $url_for_submit,
		#Content => [ 	SEQUENCE => $form_aaseq
		#		],
		Referer => $referer
		);

	# if ($create_outfile2 == 1) {
	#	print OUTFILE2 "...Printing out \$ua\n\n";
	#	print OUTFILE2 Dumper($ua);
	#	print OUTFILE2 "\n...Printing out \$request\n\n";
	#	print OUTFILE2 Dumper($request);
	# }

	my $got_timeout = 1;
	my $num_retries = 1;
	my $response;

	my $all_ok = 0;
	my $all_not_ok_msg = '';

	while (($got_timeout == 1) && ($num_retries <= $max_num_retries_of_query1)) {

		$response = $ua->request($request);

		# if ($create_outfile2 == 1) {
		#	print OUTFILE2 "\n...Printing out \$response\n\n";
		#	print OUTFILE2 Dumper($response);
		#	print OUTFILE2 "\n";
		# }

		if ($response->message eq '500') {
			$got_timeout = 1;
		} else {
			$got_timeout = 0;
		}

		if ($create_outfile2 == 1) {
			my $print_response_message = $response->message;
			print OUTFILE2 "...query 1 : print_response_message=$print_response_message, got_timeout=$got_timeout, num_retries=$num_retries\n";
		}

		$num_retries += 1;
	}

	if ($got_timeout == 0) {

		my $content = $response->content;

		# if ($create_outfile2 == 1) {
		#	print OUTFILE2 "\n...Printing out \$content\n\n";
		#	print OUTFILE2 "$content\n\n";
		# }


		# http://www.rcsb.org/pdb/explore/explore.do?structureId=1zoy

#	<span class="se_key">Method:&nbsp;&nbsp;</span>
#
#	X-RAY DIFFRACTION
#	
#	<br />
#
#	<span class="se_key">Exp. Data:</span><br />

#			<tr class="se_subitem2" >
#
#				<td width="50%">
#			
#					<span class="se_key">Resolution[&Aring;]:</span> 
#					
#					<a href='/pdb/statistics/histogram.do?mdcat=refine&amp;mditem=ls_d_res_high&amp;minLabel=0&amp;maxLabel=5&amp;numOfbars=10'>
#
#						<span class="iconSet-main icon-chart" title="View a Histogram of Resolutions">&nbsp;</span>
#
#					</a>
#
#				</td>
#				
#				<td id="se_xrayResolution">
#				
#					2.40
#				
#				</td>
#				
#			</tr>

		# http://www.rcsb.org/pdb/explore/explore.do?structureId=2LBV

#	<div id="se_unitExperimentalMethod"   >	
#	
#
#
#	<div class="lb_body3">
#
#		
#			
#			<div class="box_contentWrapper">
#
#		
#
#
#			
#	<span class="se_key">Method:&nbsp;&nbsp;</span>
#
#	SOLUTION NMR
#	
#	<br />
#
#	<span class="se_key">Exp. Data:</span><br />



		my @lines = split(/\n/, $content);

		my $method = '';
		my $resolution = '';

		my $i = 0;
		foreach my $line (@lines) {
			chomp $line;
			$line = trim($line);

			#print "$line\n";

			if ($line =~ /<span class="se_key">Method:&nbsp;&nbsp;<\/span>/) {
				my $j = $i + 2;
				$method = trim( $lines[$j] );
			}

			if ($line =~ /<td id="se_xrayResolution">/) { # X-RAY DIFFRACTION
				my $j = $i + 2;
				$resolution = trim( $lines[$j] );
			}

			if ($line =~ /<span class="se_key">EM Resolution/) { # ELECTRON MICROSCOPY, ELECTRON CRYSTALLOGRAPHY
				my $j = $i + 2;
				$resolution = trim( $lines[$j] );
			}

			# no resolution : FIBER DIFFRACTION, SOLID-STATE NMR, SOLUTION NMR

			$i++;
		}

		if ($method ne '') {
			my $output_line = $pdbid_and_chain . ':::::METHOD,,,,,' . $method . ';;;;;';
			if ($resolution ne '') {
				$output_line .= 'RESOLUTION,,,,,' . $resolution . ';;;;;';
			}
			$output_line .= ':::::' . "\n";
			print OUTFILE $output_line;
		}

	} else {
		print "     Had $max_num_retries_of_query1 timeouts - got no results.\n";
		$all_not_ok_msg = "Had $max_num_retries_of_query1 timeouts - got no results.";
	}

	if ($all_ok == 0) {
		print REDOFILE "\>$pdbid_and_chain\n\n";
	}

	close OUTFILE;
	open( OUTFILE, ">>$output_file") or
		die "Cannot open $output_file for appending\n";
}

##################################################
# Perl trim function to remove whitespace from the start and end of the string
sub trim {
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}

