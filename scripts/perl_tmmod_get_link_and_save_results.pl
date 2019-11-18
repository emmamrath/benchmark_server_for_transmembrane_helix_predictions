#!/usr/bin/perl   -w

# perl perl_tmmod_get_link_and_save_results.pl -infile alpha_polytopic_bitopic_noblanks.1line.tmmod_link

# will produce output file xxx.tmmod_results
# will produce output file xxx.tmmod_html2
# will produce output file xxx.tmmod_redo2

# The input to this program is a file containing multiple fasta sequences.
# Each fasta sequence was submitted to the MEMSAT-SVM server at http://bioinf.cs.ucl.ac.uk/psipred/?program=svmmemsat
# and the HTML results gave a link that was saved into the input file of this program.
# This program will do an HTTP GET of that link, to get the actual results.
# The REDO file looks just like the input file, and contains the sequences that always got rc 500 timeout,

# The TMpred program makes a prediction of membrane-spanning regions and their orientation. 
# The algorithm is based on the statistical analysis of TMbase, 
# a database of naturally occuring transmembrane proteins. 
# The prediction is made using a combination of several weight-matrices for scoring.



use strict;
use Getopt::Long;
use Bio::SeqIO;
use LWP;
use HTTP::Request::Common;
use Data::Dumper;

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
my $url_for_submit = 'http://liao.cis.udel.edu/website/servers/TMMOD/scripts/frame.php';

# global variables

my $max_num_retries_of_query1 = 10;
my $count_input_seqs = 0;
my $output_info;
my $referer = 'Client User Agent';

my $input_file;
GetOptions( "infile=s" => \$input_file );
if (!defined $input_file) {
	die "Usage $0 -infile 1LINEfile\n";
}

# Open output files

my $output_file = "$input_file.tmmod_results";
my $output_file_2 = "$input_file.tmmod_html2";
my $redo_file = "$input_file.tmmod_redo2";

open( OUTFILE, ">$output_file") or
	die "Cannot open $output_file for writing\n";
if ($create_outfile2 == 1) {
	open( OUTFILE2, ">$output_file_2") or
		die "Cannot open $output_file for writing\n";
}
open( REDOFILE, ">$redo_file") or
	die "Cannot open $output_file for writing\n";

# processing

open INFILE, $input_file or die $!;
my @input_lines = <INFILE>;

foreach my $input_line (@input_lines) {

	chomp( $input_line );
	$input_line = trim($input_line);
	if ($input_line ne '') {

		# write the fasta-id and fasta-sequence to the temporary file,
		# so it can be input to the tmhmm command

		my @bits = split(/:::::/, $input_line);

		my $this_id = $bits[0];
		my $this_link = $bits[1];

		if ($this_link ne '-') {

			($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
			$year += 1900;
			$mon += 1;
			$mon = sprintf("%02d", $mon);
			$mday = sprintf("%02d", $mday);
			$hour = sprintf("%02d", $hour);
			$min = sprintf("%02d", $min);
			$sec = sprintf("%02d", $sec);
			$time_str = "$year-$mon-$mday.$hour:$min:$sec";

			$count_input_seqs++;
			print         "$time_str processing input number $count_input_seqs, sequence $this_id...\n";
			# print OUTFILE "$time_str Processing $this_id\n";
			# print OUTFILE "$input_aaseq\n";
			if ($create_outfile2 == 1) {
				print OUTFILE2 "$time_str Processing $this_id\n\n";
			}

			my $ua = LWP::UserAgent->new;
			#my $newagent = 'test search program (' . $ua->agent . ')';
			my $newagent = 'Client User Agent';
			$ua->agent($newagent);

			my $request = HTTP::Request->new( GET => $this_link );
			my $response = $ua->request($request);

#<table cellpadding="0" cellspacing="
#0"><tr> <td ALIGN=RIGHT VALIGN=TOP > <font size= 2><table cellpadding="0" cellspacing="
#0"><tr> <td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Name <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 1s5h_C <font></td></tr><tr> <td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Annotation <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TM PROTEIN <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr><tr> <td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Length <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 124 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr><tr> <td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Number of predicted TMHs <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 2 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr><tr> <td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Exp number of AAs in TMHs <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 42.907185 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr><tr> <td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Exp number, first 60 AAs <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 19.532410 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr><tr> <td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Total prob of N-in <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 0.463631 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; inside <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 1 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 28 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TMhelix <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 29 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 49 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; outside <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 50 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 89 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TMhelix <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 90 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 115 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; inside <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 116 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 124 <font></td></tr><tr> <td ALIGN=LEFT VALIGN=CENTER > <font size= 2> <a href="frame.php?p=posterior&jobid=1633333927&index=0".> show posterior probabilities </a>  <font></td> </tr></table><font></td> <td ALIGN=RIGHT VALIGN=TOP > <font size= 2>&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  <font></td></tr></table><br><hr><br><br><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp;  <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp;  <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr>
#</td>

			if ($response->is_success) {

				my $output_line = $this_id . ':::::';
				print OUTFILE "$output_line\n";

				my $http_result2 = $response->content;
				my @lines2 = split(/\n/, $http_result2);

				my $need_to_output_current_result_lines = 0;

				foreach my $line2 (@lines2) {

					$line2 = trim($line2);
					if ($line2 ne '') {

						if (index($line2,"&nbsp;&nbsp;&nbsp;&nbsp;") > -1) {
							print OUTFILE "$line2\n";
						}
					}
				}

				close OUTFILE;
				open( OUTFILE, ">>$output_file") or
					die "Cannot open $output_file for writing\n";
			}
		}
	}
}

print OUTFILE "\n";
close OUTFILE;
if ($create_outfile2 == 1) {
	print OUTFILE2 "\n";
	close OUTFILE2;
}
close REDOFILE;



sub find_leading_trailing_GGG_HHH_AAA {

	# remove leading AAA, remove leading GGG, remove trailing HHH

	my $in_aaseq = shift;
	my %seq_and_index_array;
	my $seq_and_index = \%seq_and_index_array;
	$seq_and_index->{'seq'} = '';
	$seq_and_index->{'index'} = 0;
	my $aaseq = uc( $in_aaseq );

	my @aaseq_array = split(//, $aaseq);
	my $in_AAA = 1;
	my $start_of_non_AAA = 0;
	for ( my $i = 0; $i <@aaseq_array; $i++ ) {
		if ($in_AAA == 1) {
			if ($aaseq_array[$i] ne 'A') {
				$start_of_non_AAA = $i;
				$in_AAA = 0;
			}
		}
	}
	my $new_aaseq = $aaseq;
	if ($start_of_non_AAA > 0) {
		$new_aaseq = substr( $aaseq, $start_of_non_AAA );
	}
	$aaseq = $new_aaseq;

	@aaseq_array = split(//, $aaseq);
	my $in_GGG = 1;
	my $start_of_non_GGG = 0;
	for ( my $i = 0; $i <@aaseq_array; $i++ ) {
		if ($in_GGG == 1) {
			if ($aaseq_array[$i] ne 'G') {
				$start_of_non_GGG = $i;
				$in_GGG = 0;
			}
		}
	}
	$new_aaseq = $aaseq;
	if ($start_of_non_GGG > 0) {
		$new_aaseq = substr( $aaseq, $start_of_non_GGG );
	}
	$aaseq = $new_aaseq;

	@aaseq_array = split(//, $aaseq);
	my $in_HHH = 1;
	my $start_of_HHH = @aaseq_array;
	for ( my $i = @aaseq_array - 1; $i >= 0; $i-- ) {
		if ($in_HHH == 1) {
			if ($aaseq_array[$i] eq 'H') {
				$start_of_HHH = $i;
			} else {
				$in_HHH = 0;
			}
		}
	}
	$new_aaseq = $aaseq;
	if ($start_of_HHH < @aaseq_array) {
		my $allow_num_consecutive_H = 1; # allow this number of H at the end without removing them. if any more, then remove all tail HHH.
		if ($start_of_HHH <= (@aaseq_array - ($allow_num_consecutive_H + 1)) ) {
			$new_aaseq = substr( $aaseq, 0, $start_of_HHH );
		}
	}
	$aaseq = $new_aaseq;

	$seq_and_index->{'seq'} = $aaseq;
	$seq_and_index->{'index'} = index( $in_aaseq, $aaseq );
	return $seq_and_index;
}



sub trim {
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}

