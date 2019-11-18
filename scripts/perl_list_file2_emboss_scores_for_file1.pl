#!/usr/bin/perl   -w

# perl ../perl_list_file2_emboss_scores_for_file1.pl -infile list11_not_in_mpstruc_pass1_pass2_pass3.pdb_seq -align needle -infile2 one2_mpstruc_category.id.pdb_seq.needle.similarity_scores.30.binary_matrix.nonhomologous_list.pdb_seq -outfile list11_not_in_mpstruc_pass1_pass2_pass3.needle.235_mpstruc_ids


# Example of input file:
# 1a11_A:::::GSEKMSTAISVLLAQAVFLLLTSQR:::::
# 1a91_A:::::MENLNMDLLYMAAAVMMGLAAIGAAIGIGILGGKFLEGAARQPDLIPLLRTQFFIVMGLVDAIPMIAVGLGLYVMFAVA:::::
# 1afo_A:::::VQLAHHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIKK:::::

# Example of output file:
# - : 1a11_A,1a91_A,1afo_A
# 1h68_A : 28.5,53.3,47.4,
# 2ei4_A : 45.7,43.2,71.1,
# 1py6_A : 34.4,47.6,54.4,

use warnings;
use strict;
use diagnostics;
use Getopt::Long;
# use Data::Dumper;

# Check to see if -infile was specified
my ($input_file, $input_align, $input_file_2, $output_file);
GetOptions( "infile=s" => \$input_file, "align=s" => \$input_align, "infile2=s" => \$input_file_2,  "outfile=s" => \$output_file );
if ((!defined $input_file) || (!defined $input_align) || (!defined $input_file_2) || (!defined $output_file)) {
	die "Usage $0 -infile INPUTFILE\n";
}
$input_align = lc($input_align);

my $command_input_file_1 = "scrap2_temporary_in_1_$input_align.txt";
my $command_input_file_2 = "scrap2_temporary_in_2_$input_align.txt";
my $command_output_file = "scrap2_temporary_out_$input_align.txt";
my $command_input_params_file = "scrap2_temporary_params_in_$input_align.txt";

# Open output files
my $output_file_similarity = "$output_file.similarity_scores";
open( OUTFILE_SIM, ">$output_file_similarity") or
	die "Cannot open $output_file_similarity for writing\n";
my $output_file_identity = "$output_file.identity_scores";
open( OUTFILE_IDENT, ">$output_file_identity") or
	die "Cannot open $output_file_identity for writing\n";

open( PARAMS, ">$command_input_params_file") or
	die "Cannot open $command_input_params_file for writing\n";
my $params_lines = "\n\n";
print PARAMS $params_lines;
close PARAMS;

open INFILE, $input_file or die $!;
my @input_lines = <INFILE>;

open INFILE2, $input_file_2 or die $!;
my @input_lines_2 = <INFILE2>;

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

my @input_seqid_2;
my @input_aaseq_2;
foreach my $input_line_2 (@input_lines_2) {

	chomp($input_line_2);
	if ($input_line_2 ne '') {

		my @bits = split(/:::::/, $input_line_2);
		my $seqid = $bits[0];
		my $aaseq = $bits[1];
		push( @input_seqid_2, $seqid );
		push( @input_aaseq_2, $aaseq );
	}
}

my $output_line_similarity = '- : ';
my $output_line_identity = '- : ';
for ( my $i = 0; $i < @input_seqid_2; $i++ ) {
	$output_line_similarity .= $input_seqid_2[$i] . ',';
	$output_line_identity .= $input_seqid_2[$i] . ',';
}
print OUTFILE_SIM "$output_line_similarity\n";
close OUTFILE_SIM;
open( OUTFILE_SIM, ">>$output_file_similarity") or
	die "Cannot open $output_file_similarity for writing\n";
print OUTFILE_IDENT "$output_line_identity\n";
close OUTFILE_IDENT;
open( OUTFILE_IDENT, ">>$output_file_identity") or
	die "Cannot open $output_file_identity for writing\n";

for ( my $i = 0; $i < @input_seqid; $i++ ) {
	$output_line_similarity = $input_seqid[$i] . ' : ';
	$output_line_identity = $input_seqid[$i] . ' : ';
	for ( my $j = ($i + 1); $j < @input_seqid_2; $j++ ) {

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
		my $seq2 = '>' . $input_seqid_2[$j] . "\n";
		$seq2 .= $input_aaseq_2[$j] . "\n";
		print SEQ2 $seq2;
		close SEQ2;

		print "\n\nprocessing $input_align for $input_seqid[$i] and $input_seqid_2[$j]...\n";

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


