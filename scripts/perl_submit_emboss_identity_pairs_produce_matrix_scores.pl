#!/usr/bin/perl   -w

# perl perl_submit_emboss_identity_pairs_produce_matrix_scores.pl -infile opm_feb2012.unique_seq_opm_in_pdb -align water
# perl perl_submit_emboss_identity_pairs_produce_matrix_scores.pl -infile opm_feb2012.unique_seq_opm_in_pdb -align needle

# will produce output file xxx.similarity_scores and identity.identity_scores
# eg. opm_feb2012.unique_seq_opm_in_pdb.similarity_scores and opm_feb2012.unique_seq_opm_in_pdb.identity_scores

# Will output a matrix of whether each pair in a list is similar or not similar.
# Similar or not is decided by emboss water, less than 30% similar.
# Matrix output :
#	_	ignore
#	1	not similar
#	0	similar

# Example of input file:
# 1a11_A:::::GSEKMSTAISVLLAQAVFLLLTSQR:::::
# 1a91_A:::::MENLNMDLLYMAAAVMMGLAAIGAAIGIGILGGKFLEGAARQPDLIPLLRTQFFIVMGLVDAIPMIAVGLGLYVMFAVA:::::
# 1afo_A:::::VQLAHHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIKK:::::

# Example of output file:
# 1h68_A : ,53.3,47.4,
# 2ei4_A : ,,71.1,
# 1py6_A : ,,,

use warnings;
use strict;
use diagnostics;
use Getopt::Long;
# use Data::Dumper;

# Check to see if -infile was specified
my ($input_file, $input_align);
GetOptions( "infile=s" => \$input_file, "align=s" => \$input_align );
if ((!defined $input_file) || (!defined $input_align)) {
	die "Usage $0 -infile INPUTFILE\n";
}
$input_align = lc($input_align);

my $command_input_file_1 = "scrap_temporary_in_1_$input_align.txt";
my $command_input_file_2 = "scrap_temporary_in_2_$input_align.txt";
my $command_output_file = "scrap_temporary_out_$input_align.txt";
my $command_input_params_file = "scrap_temporary_params_in_$input_align.txt";

# Open output files
my $output_file_similarity = "$input_file.$input_align.similarity_scores";
open( OUTFILE_SIM, ">$output_file_similarity") or
	die "Cannot open $output_file_similarity for writing\n";
my $output_file_identity = "$input_file.$input_align.identity_scores";
open( OUTFILE_IDENT, ">$output_file_identity") or
	die "Cannot open $output_file_identity for writing\n";

open( PARAMS, ">$command_input_params_file") or
	die "Cannot open $command_input_params_file for writing\n";
my $params_lines = "\n\n";
print PARAMS $params_lines;
close PARAMS;

open INFILE, $input_file or die $!;
my @input_lines = <INFILE>;

my @input_seqid;
my @input_aaseq;
foreach my $input_line (@input_lines) {

	chomp($input_line);
	if ($input_line ne '') {

		my @bits = split(/:::::/, $input_line);
		my $seqid = $bits[0];
		my $aaseq = $bits[1];
		push( @input_seqid, $seqid );
		push( @input_aaseq, $aaseq );
	}
}

for ( my $i = 0; $i < @input_seqid; $i++ ) {
	my $output_line_similarity = $input_seqid[$i] . ' : ';
	my $output_line_identity = $input_seqid[$i] . ' : ';
	for ( my $j = 0; $j <= $i; $j++ ) {
		$output_line_similarity .= ',';
		$output_line_identity .= ',';
	}
	for ( my $j = ($i + 1); $j < @input_seqid; $j++ ) {

		#=======================================
		#
		# Aligned_sequences: 2
		# 1: scrap1
		# 2: scrap2
		# Matrix: EBLOSUM62
		# Gap_penalty: 10.0
		# Extend_penalty: 0.5
		#
		# Length: 48
		# Identity:      33/48 (68.8%)
		# Similarity:    39/48 (81.2%)
		# Gaps:           5/48 (10.4%)
		# Score: 143.5
		# 
		#
		#=======================================

		open( SEQ1, ">$command_input_file_1") or
			die "Cannot open $command_input_file_1 for writing\n";
		my $seq1 = '>' . $input_seqid[$i] . "\n";
		$seq1 .= $input_aaseq[$i] . "\n";
		print SEQ1 $seq1;
		close SEQ1;

		open( SEQ2, ">$command_input_file_2") or
			die "Cannot open $command_input_file_2 for writing\n";
		my $seq2 = '>' . $input_seqid[$j] . "\n";
		$seq2 .= $input_aaseq[$j] . "\n";
		print SEQ2 $seq2;
		close SEQ2;

		print "\n\nprocessing $input_align for $input_seqid[$i] and $input_seqid[$j]...\n";

		my $command = '/home/emma/Emmas_files_on_linuxbox/software_to_install/emboss/emboss/emboss/' . $input_align . ' -asequence ' . $command_input_file_1 . ' -bsequence ' . $command_input_file_2 . ' -outfile ' . $command_output_file . ' < ' . $command_input_params_file;
		my $command_output = `$command`; 
		#print "command = $command\n";

		open COMFILE, $command_output_file or die $!;
		my @input_lines = <COMFILE>;

		# my $num_output_lines = @command_output;
		#print "num_output_lines = $num_output_lines\n";

		my $got_all_needed_info = 0;
		my $alignment_identity = -1;
		my $alignment_similarity = -1;
		my $alignment_gaps = -1;
		foreach my $input_line (@input_lines) {
			if ($got_all_needed_info == 0) {
				chomp( $input_line );
				if ($input_line ne "") {
					my $key_to_this_line = substr($input_line, 0, 13);

					if ($key_to_this_line eq '# Identity:  ') {
						my @bits = split(/\(/, $input_line);
						my @bits2 = split(/\%/, $bits[1]);
						$alignment_identity = $bits2[0];
					} elsif ($key_to_this_line eq '# Similarity:') {
						my @bits = split(/\(/, $input_line);
						my @bits2 = split(/\%/, $bits[1]);
						$alignment_similarity = $bits2[0];
					} elsif ($key_to_this_line eq '# Gaps:      ') {
						my @bits = split(/\(/, $input_line);
						my @bits2 = split(/\%/, $bits[1]);
						$alignment_gaps = $bits2[0];
					}

					if (($alignment_identity != -1) && ($alignment_similarity != -1) && ($alignment_gaps != -1)) {
						$got_all_needed_info = 1;
					}
				}
			}
		}

		$output_line_similarity .= $alignment_similarity . ',';
		$output_line_identity .= $alignment_identity . ',';
	}
	print OUTFILE_SIM "$output_line_similarity\n";
	close OUTFILE_SIM;
	open( OUTFILE_SIM, ">>$output_file_similarity") or
		die "Cannot open $output_file_similarity for writing\n";
	print OUTFILE_IDENT "$output_line_identity\n";
	close OUTFILE_IDENT;
	open( OUTFILE_IDENT, ">>$output_file_identity") or
		die "Cannot open $output_file_identity for writing\n";
	print "\n";
}
close OUTFILE_SIM;
close OUTFILE_IDENT;


