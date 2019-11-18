#!/usr/bin/perl   -w

# perl perl_toppred_fasta_submit_local.pl -infile scrap.fasta
# perl perl_toppred_fasta_submit_local.pl -infile scrap.fasta -e
# perl perl_toppred_fasta_submit_local.pl -infile scrap.fasta -e -rmvggg

# If the -rmvggg flag is present, then remove leading and trailing GGG, HHH and AAA before submitting to prediction tool.

# If the -e flag is not present, then the TOPPRED procaryote default will be used.
# If the -e flag is present, then the TOPPRED eucaryote default will be used.

# will produce output file xxx.toppred # will also produce a temporary file xxx.toppred_temp that will be left there permanently (eg. scrap1_seqids.unique.1line.tmhmm_temp)
# will also produce a temporary file xxx.toppred_temp2 that will be left there permanently (eg. scrap1_seqids.unique.1line.tmhmm_temp2)

# The input to this program is a file containing seqids, fasta-ids, and aa-seqs, 1 line per seqid.
# The output has, for each seqid, the toppred results.
# The call to toppred needs a input file containing the 1 fasta sequence,
# which is why this program recreates the temporary file for each sequence processed.
# The second temporary file will contain the output returned by tmhmm.

# Example of input file:
# AAA26845:::::>gi|153550|gb|AAA26845.1| penicillin-binding protein (PBPB2):::::AVIASISKEMPGISISTSWDRKILETSLSSIVGSVSSEKAGLPAEEAETYLKKGYSLNDRVGTSYLEKQYEETLQGKRSVKEIHLDKYGNMESVENIEDGTKGNNIKLTIDLSFQDSVDALLKSYFNSELGNGGAKYSEGVYAVALNPKTGAVLSMSGIKHDLKTGELTPDSLGTVTNVFVPGSVVKAATISSGWENGVLSGNQTLTDQSIVFQGSAPINSWYTAFSRPMPITAVQALEYSSNAYMVQTALGLMGQTYQPNMFVGTSNLESAMGKLRSTFGEYGLGSATGIDLPDESTGFIPKEYSFANFITNAFGQFDNYTPMQLAQYVATIANDGVRVAPRIVEGIYGNNDKGGLGGLIQQLQPTEMNKVNISDSDMSVLHQGFYQVAHGTSGLTTGRAFSNGAAVSISGKTGTAESYVAGGQEANNTNAVAYAPSDNPQIAVAVVFPHNTNLTNGVGPSIARDIINLYNQHHPMN
# AAA26872:::::>gi|153613|gb|AAA26872.1| DpnII DNA methylase:::::MKIKEIKKVTLQPFTKWTGGKRQLLPVIRELIPKTYNRYFEPFVGGGALFFDLAPKDAVINDFNAELINCYQQIKDNPQELIEILKVHQEYNSKEYYLDLRSADRDERIDMMSEVQRAARILYMLRVNFNGLYRVNSKNQFNVPYGRYKNPKIVDEELISAISVYINNNQLEIKVGDFEKAIVDVRTGDFVYFDPPYIPLSETSAFTSYTHEGFSFADQVRLRDAFKRLSDTGAYVMLSNSSSALVEELYKDFNIHYVEATRTNGAKSSSRGKISEIIVTNYEK

# Example of TOPPRED output:
# Found: 7 segments
# 
# Candidate membrane-spanning segments:
# 
#  Helix Begin - End   Score Certainity
#      1    16 - 36    2.107 Certain
#      2    53 - 73    1.874 Certain
#      3    89 - 109   0.889 Putative
#      4   114 - 134   1.441 Certain
#      5   142 - 162   2.237 Certain
#      6   179 - 199   1.960 Certain
#      7   210 - 230   0.946 Putative
# 
# Total of 4 structures are to be tested



use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use IO::File;
use XML::Parser;
use Data::Dumper;

# Check to see if -infile was specified

my $input_file;
my $input_e;
my $input_rmvggg;
GetOptions( "infile=s" => \$input_file, "e=s" => \$input_e, "rmvggg" => \$input_rmvggg );
if (!defined $input_file) {
	die "Usage $0 -infile 1LINESfile\n";
}

# Open output file
my $output_file = "$input_file.toppred";
my $temp_file = "$input_file.toppred_temp";
my $temp_file_2 = "$input_file.toppred_temp2";
open( OUTFILE, ">$output_file") or
	die "Cannot open $output_file for writing\n";

# processing

open INFILE, $input_file or die $!;
my @input_lines = <INFILE>;

my $upto = 0;

foreach my $input_line (@input_lines) {

	chomp( $input_line );
	if ($input_line ne '') {

		# write the fasta-id and fasta-sequence to the temporary file,
		# so it can be input to the tmhmm command

		my @bits = split(/:::::/, $input_line);

		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		$year += 1900;
		$mon += 1;
		$mon = sprintf("%02d", $mon);
		$mday = sprintf("%02d", $mday);
		$hour = sprintf("%02d", $hour);
		$min = sprintf("%02d", $min);
		$sec = sprintf("%02d", $sec);
		my $time_str = "$year-$mon-$mday.$hour:$min:$sec";

		$upto++;
		my $input_seqid = $bits[0];
		print "$time_str : Processing sequence number $upto : $input_seqid\n";

		my $input_fastaid = $bits[1];
		my $output_redo = '-';

		my $input_aaseq_multiline = '';
		my $max_str_len = 70;
		my $input_aaseq = $bits[2];

		# if an amino acid is marked as X, TOPPRED will produce the following message:
		# 	Warning: X: not a valid aa, skipped
		# and will not count that residue, which will mean all the predicted residue positions will be wrong
		$input_aaseq =~ s/X/G/g;

		my %seq_and_index_array;
		my $seq_and_index = \%seq_and_index_array;
		$seq_and_index->{'seq'} = $input_aaseq;
		$seq_and_index->{'index'} = 0;
		if (defined($input_rmvggg)) {
			$seq_and_index = find_leading_trailing_GGG_HHH_AAA( $input_aaseq );
		}

		my $remaining_seq = $seq_and_index->{'seq'};
		my $keep_looping = 1;
		while ($keep_looping == 1) {

			my $str_len = length $remaining_seq;
			if ($str_len <= $max_str_len) {
				$keep_looping = 0;
				$input_aaseq_multiline .= "\n" . $remaining_seq;
			} else {
				$input_aaseq_multiline .= "\n" . substr($remaining_seq, 0, $max_str_len);
				$remaining_seq = substr($remaining_seq, $max_str_len) ;
			}
		}

		my $output_lines = $input_fastaid . $input_aaseq_multiline . "\n";

		open( TEMPFILE, ">$temp_file") or
			die "Cannot open $temp_file for writing\n";
		print TEMPFILE $output_lines;
		close TEMPFILE;

		# Kingdom: procaryote
		# /usr/local/bin/toppred seq.fasta -o seq.out

		# Kingdom: eucaryote
		# /usr/local/bin/toppred seq.fasta -o seq_e.out -e

		# now run the toppred command, and capture the output in this program

		my $toppred_executable = '/usr/local/bin/toppred';
		my $toppred_command = $toppred_executable . ' ' . $temp_file . ' -o ' . $temp_file_2 . ' -g none';
		if (defined($input_e)) {
			$toppred_command .= ' -e';
		}
		my $toppred_output = `$toppred_command`;

		# now process the toppred output

		open TOPPRED, $temp_file_2 or die $!;
		my @lines = <TOPPRED>;

		my $output_info = '';

		my $in_tms_area = 0;
		my $finished_tms_area = 0;

		foreach my $line (@lines) {
			$line =~ s/\n//g;
			$line =~ s/\r//g;
			if (length($line) > 0) {

				if ($finished_tms_area == 0) {
					if ($in_tms_area == 0) {
						if (length($line) >= 18) {
							if (substr($line, 0, 18) eq ' Helix Begin - End') {
								$in_tms_area = 1;
							}
						}
					} else { # $in_tms_area == 1
						my $bit = '';
						if (length($line) >= 8) {
							$bit = substr($line, 0, 8);
						}
						if ($bit eq 'Total of') {
							$in_tms_area = 0;
							$finished_tms_area = 1;
						} else {
							$bit = uc(substr($line, 0, 4));
							my @bits = split(/[\s]+/, $line);
							my $start = $bits[2] + $seq_and_index->{'index'};
							my $stop = $bits[4] + $seq_and_index->{'index'};
							$output_info .= $start . '-' . $stop . ',';
						}
					}
				}
			}
		}

		if ($output_info eq '') {
			$output_info = '-';
		}

		my $output_line = $input_seqid . ':::::' . $input_fastaid . ':::::' . $input_aaseq . ':::::' . $output_info . ':::::';

		print OUTFILE "$output_line\n";

		#my $remove_file_bit = substr($input_fastaid,1);
		#my $remove_files = $remove_file_bit . '-*.png ' . $remove_file_bit . '.hydro';
		#my $remove_command = 'rm ' . $remove_files;
		#my $remove_output = `$remove_command`;
	}
}
close OUTFILE;


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

