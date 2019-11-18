#!/usr/bin/perl   -w

# perl ../perl_report_emboss_needle_uneven_start.pl -infile1 alpha_polytopic_bitopic_peptide.opm_seq_not_equal_pdb_seq.opm_seq -infile2 alpha_polytopic_bitopic_peptide.opm_seq_not_equal_pdb_seq.pdb_seq -printfull N -outfile alpha_polytopic_bitopic_peptide.opm_seq_not_equal_pdb_seq.unequal_start_posn
# perl ../perl_report_emboss_needle_uneven_start.pl -infile1 alpha_polytopic_bitopic_peptide.opm_seq.sorted -infile2 alpha_polytopic_bitopic_peptide.renum_pdb_seq -printfull Y -outfile alpha_polytopic_bitopic_peptide.opm_seq_not_equal_renum_pdb_seq.unequal_start_posn.printfull

# perl perl_report_emboss_needle_uneven_start.pl -infile1  [input-file-1] -infile2 [input-file-2] -printfull Y -outfile [output-file]

# Example of input file:
# 1a11_A:::::GSEKMSTAISVLLAQAVFLLLTSQR:::::
# 1a91_A:::::MENLNMDLLYMAAAVMMGLAAIGAAIGIGILGGKFLEGAARQPDLIPLLRTQFFIVMGLVDAIPMIAVGLGLYVMFAVA:::::
# 1afo_A:::::VQLAHHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIKK:::::

use warnings;
use strict;
use diagnostics;
use Getopt::Long;
# use Data::Dumper;

# Check to see if -infile was specified
my ($input_file_1, $input_file_2, $output_file, $printfull);
GetOptions( "infile1=s" => \$input_file_1, "infile2=s" => \$input_file_2, "printfull=s" => \$printfull, "outfile=s" => \$output_file );
if ((!defined $input_file_1) || (!defined $input_file_2) || (!defined $printfull) || (!defined $output_file)) {
	die "Usage $0 -infile1 INPUTFILE -infile2 INPUTFILE -printfull Y -outfile OUTPUTFILE\n";
}
$printfull = uc($printfull);

my $command_input_file_1 = "scrap3_temporary_in_1_emboss.txt";
my $command_input_file_2 = "scrap3_temporary_in_2_emboss.txt";
my $command_output_file = "scrap3_temporary_out_emboss.txt";
my $command_input_params_file = "scrap3_temporary_params_in_needle.txt";

open( PARAMS, ">$command_input_params_file") or
	die "Cannot open $command_input_params_file for writing\n";
my $params_lines = "\n\n";
print PARAMS $params_lines;
close PARAMS;

# Open output file
open( OUTFILE, ">$output_file") or
	die "Cannot open $output_file for writing\n";

open INFILE, $input_file_1 or die $!;
my @input_lines_1 = <INFILE>;
my @input_seqid_1;
my @input_aaseq_1;

open INFILE, $input_file_2 or die $!;
my @input_lines_2 = <INFILE>;
my @input_seqid_2;
my @input_aaseq_2;

foreach my $input_line_1 (@input_lines_1) {
	chomp($input_line_1);
	if ($input_line_1 ne '') {
		my @bits = split(/:::::/, $input_line_1);
		my $seqid_1 = $bits[0];
		my $aaseq_1 = $bits[1];
		$aaseq_1 =~ s/\_/X/g;
		#$aaseq_1 =~ s/\_/G/g;
		push( @input_seqid_1, $seqid_1 );
		push( @input_aaseq_1, $aaseq_1 );
	}
}

foreach my $input_line_2 (@input_lines_2) {
	chomp($input_line_2);
	if ($input_line_2 ne '') {
		my @bits = split(/:::::/, $input_line_2);
		my $seqid_2 = $bits[0];
		my $aaseq_2 = $bits[1];
		$aaseq_2 =~ s/\_/X/g;
		push( @input_seqid_2, $seqid_2 );
		push( @input_aaseq_2, $aaseq_2 );
	}
}

for ( my $i = 0; $i < @input_seqid_1; $i++ ) {

	if ($input_seqid_1[$i] ne $input_seqid_2[$i]) {
		print $input_seqid_1[$i] . " is not equal to " . $input_seqid_2[$i] . "\r\n";
	}

	open( SEQ1, ">$command_input_file_1") or
		die "Cannot open $command_input_file_1 for writing\n";
	my $seq1 = '>' . $input_seqid_1[$i] . "\n";
	$seq1 .= $input_aaseq_1[$i] . "\n";
	print SEQ1 $seq1;
	close SEQ1;

	open( SEQ2, ">$command_input_file_2") or
		die "Cannot open $command_input_file_2 for writing\n";
	my $seq2 = '>' . $input_seqid_2[$i] . "\n";
	$seq2 .= $input_aaseq_2[$i] . "\n";
	print SEQ2 $seq2;
	close SEQ2;

	print "\n\nprocessing needle for $input_seqid_1[$i] and $input_seqid_2[$i]...\n";

	my $command = '/home/emma/Emmas_files_on_linuxbox/software_to_install/emboss/emboss/emboss/needle -asequence ' . $command_input_file_1 . ' -bsequence ' . $command_input_file_2 . ' -gapopen 20.0 -outfile ' . $command_output_file . ' < ' . $command_input_params_file;
	my $command_output = `$command`; 
	#print "command = $command\n";

	open COMFILE, $command_output_file or die $!;
	my @input_lines = <COMFILE>;

	my $got_all_needed_info = 0;
	my $alignment_identity = -1;
	my $alignment_similarity = -1;
	my $alignment_gaps = -1;
	my $output_line = '';

	my $aaseq_1_line_upto = 30;
	my $aaseq_2_line_upto = 32;
	my $aaseq_1_char_upto = 0;
	my $aaseq_2_char_upto = 0;
	my $aaseq_1_beginning = 0;
	my $aaseq_2_beginning = 0;

	while (($aaseq_1_beginning == 0) or ($aaseq_2_beginning == 0)) {

		my $output_line = '';

		if ($aaseq_1_beginning == 0) {
			my @bits = split( /\s+/, $input_lines[$aaseq_1_line_upto] );
			if (@bits != 4) {
				print $input_seqid_1[$i] . " has unexpected emboss needle line $aaseq_1_line_upto : " . $input_lines[$aaseq_1_line_upto] . "\r\n";
			}
			my $this_seqid = $bits[0];
			my $start_posn = $bits[1];
			my $this_aaseq = $bits[2];
			my $end_posn = $bits[3];
			my @bits2 = split( //, $this_aaseq );
			foreach my $bit2 (@bits2) {
				if ($aaseq_1_beginning == 0) {
					$aaseq_1_char_upto += 1;
					#if (($bit2 ne '-') and ($bit2 ne 'X')) {
					if ($bit2 ne '-') {
						$aaseq_1_beginning = $aaseq_1_char_upto;
					}
				}
			}
		}

		if ($aaseq_2_beginning == 0) {
			my @bits = split( /\s+/, $input_lines[$aaseq_2_line_upto] );
			if (@bits != 4) {
				print $input_seqid_2[$i] . " has unexpected emboss needle line $aaseq_2_line_upto : " . $input_lines[$aaseq_2_line_upto] . "\n";
			}
			my $this_seqid = $bits[0];
			my $start_posn = $bits[1];
			my $this_aaseq = $bits[2];
			my $end_posn = $bits[3];
			my @bits2 = split( //, $this_aaseq );
			foreach my $bit2 (@bits2) {
				if ($aaseq_2_beginning == 0) {
					$aaseq_2_char_upto += 1;
					#if (($bit2 ne '-') and ($bit2 ne 'X')) {
					if ($bit2 ne '-') {
						$aaseq_2_beginning = $aaseq_2_char_upto;
					}
				}
			}
		}

		$aaseq_1_line_upto += 4;
		$aaseq_2_line_upto += 4;
	}

	if ($aaseq_1_beginning != $aaseq_2_beginning) {

		$output_line = $input_seqid_1[$i] . ':::::' . $aaseq_1_beginning . ':::::' . $aaseq_2_beginning . ':::::';
		$output_line .= $input_seqid_1[$i] . " sequence beginning $aaseq_1_beginning in file1 is not equal to " . $input_seqid_2[$i] . " sequence beginning $aaseq_2_beginning in file 2.\n";
		print OUTFILE "$output_line";
		close OUTFILE;
		open( OUTFILE, ">>$output_file") or
			die "Cannot open $output_file for writing\n";

		if ($printfull eq 'Y') {
			if (defined $input_file_1) {
				print OUTFILE "\r\n";
				foreach my $input_line (@input_lines) {
					if ($got_all_needed_info == 0) {
						chomp( $input_line );
						if ($input_line ne '') {

							# opm_seq            1 GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG     50
							#                                                                       
							# pdb_seq            1 --------------------------------------------------      0
							#
							# opm_seq           51 GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGHGKSKL    100
							#                                                                  ||||||
							# pdb_seq            1 --------------------------------------------HGKSKL      6

							$output_line .= $input_line . "\r\n";
						}
					}
				}
				print OUTFILE "$output_line\r\n";
				close OUTFILE;
				open( OUTFILE, ">>$output_file") or
					die "Cannot open $output_file for writing\n";
			}
		}
	}

	if ($printfull eq 'Y') {
		if (defined $input_file_1) {
			print OUTFILE "\r\n\r\n";
			close OUTFILE;
			open( OUTFILE, ">>$output_file") or
				die "Cannot open $output_file for writing\n";
		}
	}

	print "\n";
}
close OUTFILE;



