#!/usr/bin/perl   -w

# perl perl_emboss_pepinfo_fasta_submit_local.pl -infile test1.1line -window 7 -cutoff -0.05 -method KYTE

# If the -rmvggg flag is present, then remove leading and trailing GGG, HHH and AAA before submitting to prediction tool.

# EBMOSS PEPINFO calculates the hydrophobicity plot for a sequence,
# Can calculate according to Kyte-Doolittle, OHM (Sweet & Eisenberg), or Consensus (Eisenberg et al) hydropathy parameters.
# Various recommended window sizes are 7, 11, and 19.
# Cutoff value for transmembrane prediction seems to be -5, although 0 might be better at separating 2 close segments.
# Methods are Kyte-Doolittle (KYTE), OHM (OHM), or Consensus (EISENBERG)

# will produce output file xxx.pepinfo_segments
# will also produce a temporary file xxx.pepinfo_segments_temp that will be left there permanently
# will also produce a temporary file xxx.pepinfo_segments_temp2 that will be left there permanently

# The input to this program is a file containing seqids, fasta-ids, and aa-seqs, 1 line per seqid.
# The output has, for each seqid, the transmembrane segment prediction results derived from the EMBOSS PEPINFO hydropathy values.
# The call to EMBOSS PEPINFO needs a input file containing the 1 fasta sequence,
# which is why this program recreates the temporary file for each sequence processed.
# The second temporary file will contain the output returned by EMBOSS PEPINFO.

# Example of input file:
# AAA26845:::::>gi|153550|gb|AAA26845.1| penicillin-binding protein (PBPB2):::::AVIASISKEMPGISISTSWDRKILETSLSSIVGSVSSEKAGLPAEEAETYLKKGYSLNDRVGTSYLEKQYEETLQGKRSVKEIHLDKYGNMESVENIEDGTKGNNIKLTIDLSFQDSVDALLKSYFNSELGNGGAKYSEGVYAVALNPKTGAVLSMSGIKHDLKTGELTPDSLGTVTNVFVPGSVVKAATISSGWENGVLSGNQTLTDQSIVFQGSAPINSWYTAFSRPMPITAVQALEYSSNAYMVQTALGLMGQTYQPNMFVGTSNLESAMGKLRSTFGEYGLGSATGIDLPDESTGFIPKEYSFANFITNAFGQFDNYTPMQLAQYVATIANDGVRVAPRIVEGIYGNNDKGGLGGLIQQLQPTEMNKVNISDSDMSVLHQGFYQVAHGTSGLTTGRAFSNGAAVSISGKTGTAESYVAGGQEANNTNAVAYAPSDNPQIAVAVVFPHNTNLTNGVGPSIARDIINLYNQHHPMN
# AAA26872:::::>gi|153613|gb|AAA26872.1| DpnII DNA methylase:::::MKIKEIKKVTLQPFTKWTGGKRQLLPVIRELIPKTYNRYFEPFVGGGALFFDLAPKDAVINDFNAELINCYQQIKDNPQELIEILKVHQEYNSKEYYLDLRSADRDERIDMMSEVQRAARILYMLRVNFNGLYRVNSKNQFNVPYGRYKNPKIVDEELISAISVYINNNQLEIKVGDFEKAIVDVRTGDFVYFDPPYIPLSETSAFTSYTHEGFSFADQVRLRDAFKRLSDTGAYVMLSNSSSALVEELYKDFNIHYVEATRTNGAKSSSRGKISEIIVTNYEK

# Example of EMBOSS PEPINFO output:
# Printing out Results from Kyte & Doolittle hydropathy parameters
# 
# Position  Residue			Result
#     1       T                           0.000
#     2       A                           0.000
#     3       A                           0.000
# 
# 
# 
# Printing out Results from OHM  Hydropathy parameters (Sweet & Eisenberg)
# 
# Position  Residue			Result
#     1       T                           0.000
#     2       A                           0.000
# 
# 
# 
# Printing out Results from Consensus parameters (Eisenberg et al)
# 
# Position  Residue			Result
#     1       T                           0.000
#     2       A                           0.000
#     3       A                           0.000



use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use IO::File;
use XML::Parser;
use Data::Dumper;

# Check to see if -infile was specified

my $input_file;
my $input_rmvggg;
my $input_window;
my $input_cutoff;
my $input_method;
GetOptions( "infile=s" => \$input_file, "rmvggg" => \$input_rmvggg, "window=s" => \$input_window, "cutoff=s" => \$input_cutoff, "method=s" => \$input_method );
if (!defined $input_file) {
	die "Usage $0 -infile 1LINESfile\n";
}
if (!defined $input_window) {
	$input_window = 7;
}
if (!defined $input_cutoff) {
	$input_cutoff = -0.05;
}
if (!defined $input_method) {
	$input_method = 'KYTE';
}

# Open output file
my $output_file = "$input_file.pepinfo_segments";
my $temp_file = "$input_file.pepinfo_segments_temp";
my $temp_file_2 = "$input_file.pepinfo_segments_temp2";
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

		my $input_aaseq_fix = $input_aaseq;
		$input_aaseq_fix =~ s/X/G/g; # convert invalid residue to valid, make sure program doesn't ignore residue and count residue positions incorrectly 

		my %seq_and_index_array;
		my $seq_and_index = \%seq_and_index_array;
		$seq_and_index->{'seq'} = $input_aaseq_fix;
		$seq_and_index->{'index'} = 0;
		if (defined($input_rmvggg)) {
			$seq_and_index = find_leading_trailing_GGG_HHH_AAA( $input_aaseq_fix );
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

		# now run the EMBOSS PEPINFO command, and capture the output in this program

		#my $local_executable = '/home/emma/Emmas_files_dont_delete/software_to_install/emboss/emboss/emboss/pepinfo -sequence ';
		my $local_executable = '/home/emma/Emmas_files_on_linuxbox/software_to_install/emboss/EMBOSS-6.4.0/emboss/pepinfo -sequence ';
		my $local_command = $local_executable . ' ' . $temp_file . ' -graph none -outfile ' . $temp_file_2 . ' -hwindow ' . $input_window;
		#print "$local_command\n";
		my $local_output = `$local_command`;

		# now process the EMBOSS PEPINFO output

		open INTERMEDIATE, $temp_file_2 or die $!;
		my @lines = <INTERMEDIATE>;

		my $output_info = '';

		my $method_string = 'Printing out Results from Kyte & Doolittle hydropathy parameters';
		if ($input_method eq 'OHM') {
			$method_string = 'Printing out Results from OHM  Hydropathy parameters (Sweet & Eisenberg)';
		} elsif ($input_method eq 'EISENBERG') {
			$method_string = 'Printing out Results from Consensus parameters (Eisenberg et al)';
		}
		my $heading_string = 'Position  Residue			Result';

		my $in_tms_area = 0;
		my $finished_tms_area = 0;
		my $seen_heading_string = 0;
		my $in_a_tms = 0;
		my $position = -1;

		foreach my $line (@lines) {
			$line =~ s/\n//g;
			$line =~ s/\r//g;
			my $length_line = length($line);
			if (length($line) > 0) {

				if ($finished_tms_area == 0) {
					if ($in_tms_area == 0) {
						if ($line eq $method_string) {
							$in_tms_area = 1;
						}
					} else { # $in_tms_area == 1
						if ($line eq $heading_string) {
							$seen_heading_string = 1;
						} else {
							my @bits = split(/[\s]+/, $line);
							$position = $bits[1] + $seq_and_index->{'index'};
							my $score = $bits[3];
							if ($score >= $input_cutoff) {
								if ($in_a_tms == 0) {
									$output_info .= $position . '-';
									$in_a_tms = 1;
								}
							} else {
								if ($in_a_tms == 1) {
									$output_info .= $position . ',';
									$in_a_tms = 0;
								}
							}
						}
					}
				}
			} else { # $length_line == 0
				if ($in_tms_area == 1) {
					if ($seen_heading_string == 1) {
						$in_tms_area = 0;
						$finished_tms_area = 1;
					}
				}
			}
		}
		if ($in_a_tms == 1) {
			$output_info .= $position . ',';
		}

		if ($output_info eq '') {
			$output_info = '-';
		}

		my $output_line = $input_seqid . ':::::' . $input_fastaid . ':::::' . $input_aaseq . ':::::' . $output_info . ':::::';

		print OUTFILE "$output_line\n";

		# my $remove_file_bit = substr($input_fastaid,1);
		# my $remove_files = $remove_file_bit . '-*.png ' . $remove_file_bit . '.hydro';
		# my $remove_command = 'rm ' . $remove_files;
		# my $remove_output = `$remove_command`;
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

