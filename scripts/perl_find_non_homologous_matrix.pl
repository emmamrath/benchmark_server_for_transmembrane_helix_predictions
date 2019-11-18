#!/usr/bin/perl   -w

# perl perl_find_non_homologous_matrix.pl -infile opm_transmembrane_alphahelical.some.txt.1line.noHHH.noGGG.similarity_matrix
# perl perl_find_non_homologous_matrix.pl -exclude list1.txt -infile opm_transmembrane_alphahelical.some2.similarity_matrix

# will produce output file xxx.nonhomologous_matrix (eg. opm_transmembrane_alphahelical.some.txt.1line.noHHH.nonhomologous_matrix)
# will produce output file xxx.nonhomologous_list (eg. opm_transmembrane_alphahelical.some.txt.1line.noHHH.nonhomologous_list)

# Will output a matrix of whether each pair in a list is similar or not similar.
# Similar or not is decided by emboss water, less than 30% similar.
# Matrix output :
#	X	ignore (because one or both sequences are not in the output list, or its the entry for a sequence with itself)
#	1	not similar
#	0	similar (so should not be in the output file at all)

# The algorithm used here to get a list of non-homologous sequences is algorithm 2 in :
# "Selection of representative protein data sets"
# Uwe Hobohm, Michael Scharf, Reinhard Schneider, Chris Sander
# Protein Science 1992, pp. 409-417

# If there is an exclude file, don't have those sequences in the final list.

# Example of input file:
# 1h68_A : _1101111111111111111111110101011111000111111111111001010011101001111111111111110111111111111111110111111111111111110101111111111111111111101011011110111111011111111111111111111110111111111111111110111111111011111111111110101111111111111111111111111111111111111111111010011110111011110111111111111111101111
# 2ei4_A : __011101111111111111111111011011111000111111111111001000011101111111101111011111011111111111011111111011111111111111111011111111111111011111111111111111111111111111111111111111100111111111111111011111111111111111111111111111111111111111111111111111111111111111111111101111110111111111011110111111111111101
# 1py6_A : ___11101111111111111111111111001111101111111111101111001011111111111111111011110101110111111111111110111111111111011111111111111111111111111111111111111011011111111111111111111110111111111111110111110111111111111111111110101111111111111111111111111111111111111111111111111110111111111011101111001111110111
# 1h2s_A : ____1111111111111111111110101011111000111111111111001010011101101111111111111110111111111111111110111111111111111110101111111111111111111101011011110111111011111111111111111111110111111111111111110111111111011111111111110101111111111111111111111111111111111111111111010011110111011110111111111111111101111
# 1h2s_B : _____111101111110101001110110101111000101101111101100001100010110110111011011110000001100011111100010011100100111110101111111010011110111110110011111111100111111101110111111100100111111111111110110111011111011111111111110100111111101111110111011111011111111111111110010110111100111100000010111001010101101

# Example of exclude file:
# 3dh4_A


use warnings;
use strict;
use diagnostics;
use Getopt::Long;
# use Data::Dumper;

# Check to see if -infile was specified
my $input_file;
my $exclude_file;
GetOptions( "infile=s" => \$input_file, "exclude=s" => \$exclude_file );
if (!defined $input_file) {
	die "Usage $0 -infile INPUTFILE\n";
}

# Get the list of exclusions
my @exclude_list;
if (defined $exclude_file) {
	open EXCLUDEFILE, $exclude_file or die $!;
	my @exclude_lines = <EXCLUDEFILE>;
	foreach my $input_line (@exclude_lines) {
		chomp($input_line);
		if ($input_line ne '') {
			my $exclude_item = $input_line;
			push(@exclude_list, $exclude_item);
		}
	}	
}

# Open output files
my $matrix_output_file = "$input_file.nonhomologous_matrix";
open( MATRIXOUT, ">$matrix_output_file") or
	die "Cannot open $matrix_output_file for writing\n";
my $list_output_file = "$input_file.nonhomologous_list";
open( LISTOUT, ">$list_output_file") or
	die "Cannot open $list_output_file for writing\n";

open INFILE, $input_file or die $!;
my @input_lines = <INFILE>;

# read in the input matrix that says which sequences are similar to other sequences

my @input_seqid;
my @input_matrix_line;
foreach my $input_line (@input_lines) {
	chomp($input_line);
	if ($input_line ne '') {
		my @bits = split(/ : /, $input_line);
		my $seqid = $bits[0];
		my $matrix_line = $bits[1];
		push( @input_seqid, $seqid );
		push( @input_matrix_line, $matrix_line );
	}
}

# initialise the output matrix

my @matrix;
my @orig_matrix;
for ( my $i = 0; $i < @input_seqid; $i++ ) {
	my @matrix_line_array;
	my @orig_matrix_line_array;
	for ( my $j = 0; $j < @input_seqid; $j++ ) {
		my $this_cell = '_';
		if ($i == $j) {
			# this stops the program from taking count of a sequence's similarity with itself
			# which will always be 100% and will always cause the sequence to fail similarity criteria of less than 100%
			$this_cell = 'x';
		}
		$matrix_line_array[$j] = $this_cell;
		$orig_matrix_line_array[$j] = $this_cell;
	}
	$matrix[$i] = \@matrix_line_array;
	$orig_matrix[$i] = \@orig_matrix_line_array;
}

# copy the input matrix to the output matrix, to be the starting point of the output matrix

for ( my $i = 0; $i < @input_seqid; $i++ ) {
	my @matrix_line_array = split(//, $input_matrix_line[$i]); 
	for ( my $j = 0; $j < @input_seqid; $j++ ) {
		if ($matrix_line_array[$j] eq '1') {
			$matrix[$i][$j] = '1';
			$matrix[$j][$i] = '1';
		} elsif ($matrix_line_array[$j] eq '0') {
			$matrix[$i][$j] = '0';
			$matrix[$j][$i] = '0';
		}
	}
}

# exclude any excluded sequences

foreach my $exclude_item (@exclude_list) {
	my $exclude_i = -1;
	for ( my $i = 0; $i < @input_seqid; $i++ ) {
		if ($exclude_item eq $input_seqid[$i]) {
			$exclude_i = $i;
		}
	}
	my $i = $exclude_i;
	for ( my $j = 0; $j < @input_seqid; $j++ ) {
		$matrix[$i][$j] = '0';
	}
	my $j = $exclude_i;
	for ( my $i = 0; $i < @input_seqid; $i++ ) {
		$matrix[$i][$j] = '0';
	}
}

# keep on finding and removing sequences until no more homologous sequences are left in the list

my $no_more_homologous_seqs_left = 0;

{
	# are there any pairs left in the list that still are homologous/similar?
	my $found_homologous = 0;
	for ( my $i = 0; $i < @input_seqid; $i++ ) {
		for ( my $j = 0; $j < @input_seqid; $j++ ) {
			if ($matrix[$i][$j] eq '0') {
				$found_homologous = 1;
			}
		}
	}
	if ($found_homologous == 0) {
		$no_more_homologous_seqs_left = 1;
	}
}

while ($no_more_homologous_seqs_left == 0) {

	# find and remove the sequence that has the most connections (similarity) to other sequences

	my $most_connected_i = -1;
	my $most_connected_count = -1;
	for ( my $i = 0; $i < @input_seqid; $i++ ) {
		my $this_connected_count = 0;
		for ( my $j = 0; $j < @input_seqid; $j++ ) {
			if ($matrix[$i][$j] eq '0') {
				$this_connected_count++;
			}
		}
		if ($this_connected_count > $most_connected_count) {
			$most_connected_count = $this_connected_count;
			$most_connected_i = $i;
		}
	}
	{
		my $i = $most_connected_i;
		for ( my $j = 0; $j < @input_seqid; $j++ ) {
			$matrix[$i][$j] = 'x';
		}
		my $j = $most_connected_i;
		for ( my $i = 0; $i < @input_seqid; $i++ ) {
			$matrix[$i][$j] = 'x';
		}
	}

	# are there any pairs left in the list that still are homologous/similar?

	my $found_homologous = 0;
	for ( my $i = 0; $i < @input_seqid; $i++ ) {
		for ( my $j = 0; $j < @input_seqid; $j++ ) {
			if ($matrix[$i][$j] eq '0') {
				$found_homologous = 1;
			}
		}
	}
	if ($found_homologous == 0) {
		$no_more_homologous_seqs_left = 1;
	}
}

# In a final pass, reinstate any removed sequences that don't have any neighbours in the selected set.
# In practice, this final pass rarely increases the size of the selected set.

my @not_in_list;
my @is_in_list;
for ( my $i = 0; $i < @input_seqid; $i++ ) {
	my $is_in_final_list = 0;
	for ( my $j = 0; $j < @input_seqid; $j++ ) {
		if ($matrix[$i][$j] eq '1') {
			$is_in_final_list = 1;
		}
	}
	if ($is_in_final_list == 1) {
		push(@is_in_list, $i);
	} else {
		push(@not_in_list, $i);
	}
}
for ( my $k = 0; $k < @not_in_list; $k++ ) {
	my $i = $not_in_list[$k];
	my $found_homologous = 0;
	for ( my $l = 0; $l < @is_in_list; $l++ ) {
		my $j = $is_in_list[$l];
		if ($orig_matrix[$i][$j] eq '0') {
			$found_homologous = 1;
		}
	}
	if ($found_homologous == 0) {
		push(@is_in_list, $i);
		for ( my $j = 0; $j < @input_seqid; $j++ ) {
			$matrix[$i][$j] = $orig_matrix[$i][$j];
			$matrix[$j][$i] = $orig_matrix[$j][$i];
		}
	}
}

# print output list of non-homologous sequences

for ( my $i = 0; $i < @input_seqid; $i++ ) {
	my $found_non_x = 0;
	for ( my $j = 0; $j < @input_seqid; $j++ ) {
		if ($matrix[$i][$j] eq '1') {
			$found_non_x = 1;
		}
	}
	if ($found_non_x == 1) {
		my $output_line = $input_seqid[$i];
		print LISTOUT "$output_line\n";
	}
}
close LISTOUT;

# print output matrix

for ( my $i = 0; $i < @input_seqid; $i++ ) {
	my $output_line = '';
	for ( my $j = 0; $j < @input_seqid; $j++ ) {
		$output_line .= $matrix[$i][$j];
	}
	print MATRIXOUT "$output_line\n";
}
close MATRIXOUT;



