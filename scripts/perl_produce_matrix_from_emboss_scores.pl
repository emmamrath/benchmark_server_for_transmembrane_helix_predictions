#!/usr/bin/perl   -w

# perl perl_produce_matrix_from_emboss_scores.pl -cutoff 30 -infile opm_transmembrane_alphahelical.some.txt.1line.noHHH.noGGG.noAAA.noGGG.similarity_scores
# perl perl_produce_matrix_from_emboss_scores.pl -cutoff 60 -infile scrap_a_few.txt.similarity_scores

# will produce output file xxx.cutoff.binary_matrix (eg. opm_transmembrane_alphahelical.some.txt.1line.noHHH.noGGG.noAAA.noGGG.20.similarity_scores.binary_matrix)

# Will output a matrix of whether each pair in a list is similar or not similar.
# Similar or not is decided by emboss water, less than 30% similar.
# Matrix output :
#	_	ignore
#	1	not similar
#	0	similar

# Example of input file:
# 1h68_A : ,53.3,47.4,
# 2ei4_A : ,,71.1,
# 1py6_A : ,,,

# Example of output file:
# 1h68_A : _00
# 2ei4_A : __1
# 1py6_A : ___

use warnings;
use strict;
use diagnostics;
use Getopt::Long;
# use Data::Dumper;

# Check to see if -infile was specified
my $input_file;
my $input_cutoff;
GetOptions( "infile=s" => \$input_file, "cutoff=s" => \$input_cutoff );
if (!defined $input_file) {
	die "Usage $0 -infile SCORESfile\n";
}

my $cutoff_for_alignment_similarity = 30;
if (defined $input_cutoff) {
	$cutoff_for_alignment_similarity = $input_cutoff;
}

# Open output files
my $output_file = "$input_file.$input_cutoff.binary_matrix";
open( OUTFILE, ">$output_file") or
	die "Cannot open $output_file for writing\n";

open INFILE, $input_file or die $!;
my @input_lines = <INFILE>;

my @input_seqid;
my @input_scores;
foreach my $input_line (@input_lines) {
	chomp($input_line);
	if ($input_line ne '') {
		my @bits = split(/ : /, $input_line);
		my $seqid = $bits[0];
		my $scores = $bits[1];
		push( @input_seqid, $seqid );
		push( @input_scores, $scores );
	}
}

for ( my $i = 0; $i < @input_seqid; $i++ ) {
	my $output_line = $input_seqid[$i] . ' : ';
	for ( my $j = 0; $j <= $i; $j++ ) {
		$output_line .= '_';
	}
	my @scores = split(/,/, $input_scores[$i]);
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

		my $alignment_similarity = trim($scores[$j]);

		if ($alignment_similarity >= $cutoff_for_alignment_similarity) {
			# print "\nthe proteins are similar\n";
			$output_line .= '0';
		} else {
			# print "\nthe proteins are not similar\n";
			$output_line .= '1';
		}

	}
	print OUTFILE "$output_line\n";
	close OUTFILE;
	open( OUTFILE, ">>$output_file") or
		die "Cannot open $output_file for writing\n";
}
close OUTFILE;



# Perl trim function to remove whitespace from the start and end of the string
sub trim {
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}

