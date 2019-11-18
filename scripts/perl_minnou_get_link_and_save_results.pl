#!/usr/bin/perl   -w

# perl perl_minnou_get_link_and_save_results.pl -infile test1.1line.minnou_link
# perl perl_minnou_get_link_and_save_results.pl -infile test2.1line.minnou_link
# perl perl_minnou_get_link_and_save_results.pl -infile alpha_polytopic_bitopic_noblanks.1line.minnou_link
# perl perl_minnou_get_link_and_save_results.pl -infile beta_barrel_noblanks.1line.minnou_link
# perl perl_minnou_get_link_and_save_results.pl -infile soluble.1line.minnou_link

# will produce output file xxx.minnou_results
# will produce output file xxx.minnou_html2
# will produce output file xxx.minnou_redo2

# Sequences were submitted to the minnou server http://minnou.cchmc.org/
# to predict transmembrane helices.


# Example results from server :
#
# <html>
# <head>
# <meta http-equiv="Refresh" content="0;URL=http://polyview.cchmc.org/cgi-bin/pr_picture.cgi?FName=493004&AddInfo=tm&HRN_TM=14-37,45-68,84-106,112-133,139-167,172-196,207-235">
# </head>
# </html>
# 
# <html>
# <head>
# <meta http-equiv="Refresh" content="0;URL=http://polyview.cchmc.org/cgi-bin/pr_picture.cgi?FName=493005&AddInfo=tm&HRN_TM=14-39,46-70,89-113,121-142,151-174,183-209,225-252">
# </head>
# </html>
# 
# <html>
# <head>
# <meta http-equiv="Refresh" content="0;URL=http://polyview.cchmc.org/cgi-bin/pr_picture.cgi?FName=493006&AddInfo=tm&HRN_TM=">
# </head>
# </html>



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
my $url_for_submit = '';

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

my $output_file = "$input_file.minnou_results";
my $output_file_2 = "$input_file.minnou_html2";
my $redo_file = "$input_file.minnou_redo2";

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

		my @bits = split(/:::::/, $input_line);

		my $this_id = $bits[0];
		my $this_link = $bits[1];
		my $this_seq = $bits[2];

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

# <html>
# <head>
# <meta http-equiv="Refresh" content="0;URL=http://polyview.cchmc.org/cgi-bin/pr_picture.cgi?FName=493005&AddInfo=tm&HRN_TM=14-39,46-70,89-113,121-142,151-174,183-209,225-252">
# </head>
# </html>
# 
# <html>
# <head>
# <meta http-equiv="Refresh" content="0;URL=http://polyview.cchmc.org/cgi-bin/pr_picture.cgi?FName=493006&AddInfo=tm&HRN_TM=">
# </head>
# </html>

			if ($response->is_success) {

				#my $output_line = $this_id . ':::::';
				#print OUTFILE "$output_line\n";

				my $http_result2 = $response->content;
				my @lines2 = split(/\n/, $http_result2);

				my $need_to_output_current_result_lines = 0;

				foreach my $line2 (@lines2) {

					$line2 = trim($line2);
					if ($line2 ne '') {

						if (index($line2,"HRN_TM=") > -1) {
							my @bits = split( /HRN_TM\=/, $line2 );
							my @bits2 = split( /\"\>/, $bits[1] );
							my $this_result = $bits2[0];
							if ($this_result eq '') {
								$this_result = '-';
							}
							my $outline = $this_id . ':::::' . $this_result . ':::::';
							print OUTFILE "$outline\n";
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

