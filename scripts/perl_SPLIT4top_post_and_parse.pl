#!/usr/bin/perl   -w

# perl perl_SPLIT4top_post_and_parse.pl -infile scrap.fasta -rmvggg

# If the -rmvggg flag is present, then remove leading and trailing GGG, HHH and AAA before submitting to prediction tool.

# will produce output file xxx.split4 (eg. scrap.fasta.split4)
# will produce output file xxx.split4_html (eg. scrap.fasta.split4_html)
# will produce output file xxx.split4_redo (eg. scrap.fasta.split4_redo)

# The input to this program is a file containing multiple fasta sequences.
# Each fasta sequence is submitted to the SPLIT4 server at http://split.pmfst.hr/split/4/
# The REDO file looks just like the input file, and contains the sequences that always got rc 500 timeout,

# Split 4.0 reference:
# Juretic, D., Zoranic, L., Zucic, D. "Basic charge clusters and predictions of membrane protein topology"
# J. Chem. Inf. Comput. Sci. Vol. 42, pp. 620-632, 2002. 



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
my $url_for_submit = 'http://split.pmfst.hr/cgi-bin/split/4/tamburr';

# global variables

my $max_num_retries_of_query1 = 10;
my $count_input_seqs = 0;
my $output_info;
my $referer = 'perl program';

my $fasta_file;
my $input_rmvggg;
GetOptions( "infile=s" => \$fasta_file, "rmvggg" => \$input_rmvggg );
if (!defined $fasta_file) {
	die "Usage $0 -infile FASTAfile -pat pattern\n";
}

# Open output files

my $output_file = "$fasta_file.split4";
my $output_file_2 = "$fasta_file.split4_html";
my $redo_file = "$fasta_file.split4_redo";

open( OUTFILE, ">$output_file") or
	die "Cannot open $output_file for writing\n";
if ($create_outfile2 == 1) {
	open( OUTFILE2, ">$output_file_2") or
		die "Cannot open $output_file for writing\n";
}
open( REDOFILE, ">$redo_file") or
	die "Cannot open $output_file for writing\n";

my $in  = Bio::SeqIO->new(-file => $fasta_file , '-format' => 'Fasta');

my %seq_and_index_array;
my $seq_and_index = \%seq_and_index_array;

while ( my $seq = $in->next_seq() ) {
	process_1_fasta_seq($seq);
}

print OUTFILE "\n";
close OUTFILE;
if ($create_outfile2 == 1) {
	print OUTFILE2 "\n";
	close OUTFILE2;
}
close REDOFILE;

sub process_1_fasta_seq {

	my $seq = shift;
	# print "process_1_fasta_seq: Sequence ",$seq->id," first 10 bases ",$seq->subseq(1,10),"\n";
	# print Dumper($seq);
	my $input_aaseq = $seq->seq;
	my $input_seqid = $seq->id . '';
	my $input_fastaid = '>' . $input_seqid;
	my $input_seq_desc = '';
	if (defined($seq->desc)) {
		$input_seq_desc = $seq->desc . '';
	}

	my $input_aaseq_fix = $input_aaseq;
	$input_aaseq_fix =~ s/X/G/g; # convert invalid residue to valid, make sure program doesn't ignore residue and count residue positions incorrectly 

	$seq_and_index->{'seq'} = $input_aaseq_fix;
	$seq_and_index->{'index'} = 0;
	if (defined($input_rmvggg)) {
		$seq_and_index = find_leading_trailing_GGG_HHH_AAA( $input_aaseq_fix );
	}

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
	print         "$time_str processing input number $count_input_seqs, sequence $input_seqid...\n";
	# print OUTFILE "$time_str Processing $input_seqid\n";
	# print OUTFILE "$input_aaseq\n";
	if ($create_outfile2 == 1) {
		print OUTFILE2 "$time_str Processing $input_seqid\n\n";
	}

	my $ua = LWP::UserAgent->new;
	my $newagent = 'test search program (' . $ua->agent . ')';
	$ua->agent($newagent);

	my $form_aaseq = $seq_and_index->{'seq'};

	my $request = POST ( $url_for_submit,
		Content => [ 	SEQUENCE => $form_aaseq
				],
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

		# http://split.pmfst.hr/cgi-bin/split/4/tamburr

		# <H1 ALIGN=CENTER><FONT SIZE=5>
		# <A HREF="http://split.pmfst.hr/split/4/results/results755.html">
		# Click here to see the prediction results!<BR>
		# </A>
		# </FONT></H1>

		# http://split.pmfst.hr/split/4/results/results755.html

		# <FONT SIZE=4>
		# <A HREF="../out/out755.html#">
		# Click here for numeric data (the original output of the prediction program)!
		# </A>
		# </FONT>

		# http://split.pmfst.hr/split/4/out/out755.html#

		#<PRE>
		#
		#  PREDICTED TM HELICES + TOPOLOGY PREDICTION
		#
		#  scale  #TMS   pred. topo.
		#    4     12         INTOP
		#  TM HELICES
		#    1        27          49
		#    2        70          89
		#    3       100         121
		#    4       125         143
		#    5       164         188
		#    6       193         213
		#    7       258         273                  SHORT
		#                                             AMPHI
		#    8       297         317
		#    9       327         347
		#   10       354         381
		#   11       387         412
		#   12       418         440


		#  PREDICTED TM HELICES + TOPOLOGY PREDICTION
		#
		#  scale  #TMS   pred. topo.
		#   17     10         INTOP
		#  TM HELICES
		#    1        13          34
		#    2        49          79
		#    3        97         111                  SHORT
		#    4       115         138                  AMPHI
		#    5       144         169
		#    6       188         208
		#    7       224         239                  SHORT
		#    8       243         270                  AMPHI
		#    9       275         294
		#   10       298         322


		#  PREDICTED TM HELICES + TOPOLOGY PREDICTION
		#
		#  scale  #TMS   pred. topo.
		#   52      6         OUTOP
		#  TM HELICES
		#    1        35          61
		#    2        81         108
		#    3       112         132
		#    4       150         174
		#    5       198         226
		#    6       258         285


		my @lines = split(/\n/, $content);

		$output_info = '';
		my $results_number = '';

		foreach my $line (@lines) {

			# <A HREF="http://split.pmfst.hr/split/4/results/results755.html">

			# <A HREF="http://split4.pmfst.hr/split/4/out/out333.dat">

			if ($results_number eq '') {
				if (length($line) > 0) {
					#if (substr($line,0,54) eq '<A HREF="http://split.pmfst.hr/split/4/results/results') {
					#	my $results_number_string = substr($line,54);
					#	my @bits = split( /\.html/, $results_number_string );
					#	$results_number = $bits[0];
					#	print "     results_number = $results_number     from $line\n";
					#}

					# <A HREF="http://split4.pmfst.hr/split/4/out/out333.dat"> out.dat

					if ( index( $line, 'http://split4.pmfst.hr/split/4/out/out' ) > -1 ) {
						my @bits = split( /http\:\/\/split4\.pmfst\.hr\/split\/4\/out\/out/, $line );
						my $bit1 = $bits[1];
						my @bits2 = split( /\.dat\"\>/, $bit1 );
						$results_number = $bits2[0];
						print "     results_number = $results_number     from $line\n";
					}
				}
			}
		}

		if ($results_number ne '') {

			my $output_topology = '-';

			#my $url2 = 'http://split.pmfst.hr/split/4/out/out' . $results_number . '.html#';
			my $url2 = 'http://split4.pmfst.hr/split/4/out/out' . $results_number . '.html#';
			my $request2 = HTTP::Request->new( GET => $url2 );
			my $response2 = $ua->request($request2);

			if ($response2->is_success) {
				my $http_result2 = $response2->content;
				my @lines2 = split(/\n/, $http_result2);

				my $found_output_in_content = 0;
				my $seen_tm_helices_line = 0;
				my $seen_motif_bias_line = 0;

				my $i = 0;
				foreach my $line (@lines2) {
					#print "$line\n";#####

					if (length($line) >= 44) {
						#if (substr($line,0,44) eq '  PREDICTED TM HELICES + TOPOLOGY PREDICTION') {
						if ( index( $line, 'PREDICTED TM HELICES + TOPOLOGY PREDICTION' ) > -1 ) {
							$found_output_in_content = 1;
							my $j = $i + 3;
							my $line3 = $lines2[$j];
							if ( index( $line3, 'INTOP' ) > -1 ) {
								$output_topology = 'i';
							} elsif ( index( $line3, 'OUTOP' ) > -1 ) {
								$output_topology = 'o';
							}
							#print "output_topology=$output_topology\n";#####
						}
					}

					if (length($line) >= 12) {
						#if (substr($line,0,12) eq '  motif bias') {
						if ( index( $line, 'motif bias' ) > -1 ) {
							$seen_motif_bias_line = 1;
						}
						#print "seen_motif_bias_line=$seen_motif_bias_line\n";#####
					}

					if (($seen_tm_helices_line == 1) && ($seen_motif_bias_line == 0)) {
						$line = trim($line);
						if (length($line) > 0) {
							my $ignore_this_line = 0;
							$line =~ s/(\s)+/ /g;
							if ( index( $line, 'SHORT' ) > -1 ) {
								$ignore_this_line = 1;
							} elsif ( index( $line, 'AMPHI' ) > -1 ) {
								$ignore_this_line = 1;
							}
							my @bits2 = split( /(\s)+/, $line );
							if ($#bits2 < 4) {
								$ignore_this_line = 1;
							}
							if ($ignore_this_line == 0) {
								my $from = $bits2[2] + $seq_and_index->{'index'};
								my $to = $bits2[4] + $seq_and_index->{'index'};
								$output_info .= $from . '-' . $to . ',';
								#print "output_info=$output_info\n";#####
							}
						}
					}

					if (length($line) >= 12) {
						#if (substr($line,0,12) eq '  TM HELICES') {
						if (( index( $line, 'TM HELICES' ) > -1 ) && ( index( $line, 'PREDICTED TM HELICES + TOPOLOGY PREDICTION' ) == -1 )) {
							$seen_tm_helices_line = 1;
							#print "seen_tm_helices_line=$seen_tm_helices_line\n";#####
						}
					}
					$i++;
				}

				if ($found_output_in_content == 1) {
					#print "all_ok=$all_ok\n";#####
					$all_ok = 1;
				}

				if ($all_ok == 1) {

					if ($output_info eq '') {
						$output_info = '-';
					}

					my $output_line = $input_seqid . ':::::' . $input_aaseq . ':::::' . $output_info . ':::::' . $output_topology . ':::::';

					print OUTFILE "$output_line\n";
					#print "$output_line\n";#####
					close OUTFILE;
					open( OUTFILE, ">>$output_file") or
						die "Cannot open $output_file for writing\n";

				} else {
					print "     Couldn't find results in the output from the second submit.\n";
					$all_not_ok_msg = "Couldn't find results in the output from the second submit.";
				}

			} else {
				print "     Didn't get output from the second submit.\n";
				$all_not_ok_msg = "Didn't get output from the second submit.";
			}

		} else {
			print "     Didn't get output from the first submit.\n";
			$all_not_ok_msg = "Didn't get output from the first submit.";
		}

	} else {
		print "     Had $max_num_retries_of_query1 timeouts - got no results.\n";
		$all_not_ok_msg = "Had $max_num_retries_of_query1 timeouts - got no results.";
	}

	if ($all_ok == 0) {
		print REDOFILE "\>$input_seqid\n$input_aaseq\n\n";
	}
}

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


# Perl trim function to remove whitespace from the start and end of the string
sub trim {
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}



##### This program doesn't do the following yet : 

# 1bcp_C :

#  PREDICTED TM HELICES + TOPOLOGY PREDICTION
#
#  scale  #TMS   pred. topo.
#   52      1         INTOP
#  TM HELICES
#    1       149         163                  SHORT
#                                             AMPHI
#
#  motif bias:          0
#  charge bias:         2
#  charge difference:   4
#  overall bias:        6


# 1pw4_A :

#  PREDICTED TM HELICES + TOPOLOGY PREDICTION
#
#  scale  #TMS   pred. topo.
#    4     12         INTOP
#  TM HELICES
#    1        27          49
#    2        70          89
#    3       100         121
#    4       125         143
#    5       164         188
#    6       193         213
#    7       258         273                  SHORT
#                                             AMPHI
#    8       297         317
#    9       327         347
#   10       354         381
#   11       387         412
#   12       418         440
#
#  motif bias:          4
#  charge bias:        16
#  charge difference:   1
#  overall bias:       29


