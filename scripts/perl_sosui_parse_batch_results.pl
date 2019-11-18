#!/usr/bin/perl   -w

# perl perl_sosui_parse_batch_results.pl -infile identity_30_graded.1line -infile2 identity_30_graded.sosui_batch_results.txt

# will produce output file xxx.sosui (eg. identity_30_graded.1line.sosui)

# Sequences were submitted to the SOSUI batch server http://bp.nuap.nagoya-u.ac.jp/sosui/sosuiG/sosuigsubmit.html
# to predict transmembrane helices.
# This program reads the sequences and reads the batch output,
# and produces a file where each line contains the sequence and the predictions.

# Example of infile :
# 1uaz_A:::::>1uaz_A:::::TAAVGADLLGDGRPETLWLGIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIGLTEVQVGSEMLDIYYARYADWLFTTPLLLLDLALLAKVDRVSIGTLVGVDALMIVTGLVGALSHTPLARYTWWLFSTICMIVVLYFLATSLRAAAKERGPEVASTFNTLTALVLVLWTAYPILWIIGTEGAGVVGLGIETLLFMVLDVTAKVGFGFILLRSRAILGDTEAPEPSAGAEASAAD:::::
# 3ddl_A:::::>3ddl_A:::::MLQELPTLTPGQYSLVFNMFSFTVATMTASFVFFVLARNNVAPKYRISMMVSALVVFIAGYHYFRITSSWEAAYALQNGMYQPTGELFNDAYRYVDWLLTVPLLTVELVLVMGLPKNERGPLAAKLGFLAALMIVLGYPGEVSENAALFGTRGLWGFLSTIPFVWILYILFTQLGDTIQRQSSRVSTLLGNARLLLLATWGFYPIAYMIPMAFPEAFPSNTPGTIVALQVGYTIADVLAKAGYGVLIYNIAKAKSEEEGFNVSEMVEPATASA:::::

# Example of infile2 :
# >1uaz_A 	MEMBRANE PROTEIN 	7
# REGION 17-39,50-72,89-111,113-135,142-164,177-199,209-231
# >3ddl_A 	MEMBRANE PROTEIN 	7
# REGION 13-35,46-65,93-114,125-147,155-177,192-214,226-248
# >1bcc_N 	SOLUBLE PROTEIN 	NONE
# >2fyn_B 	MEMBRANE PROTEIN 	1
# REGION 228-248

# References :
#
# Hirokawa T., Boon-Chieng S., and Mitaku S., Bioinformatics, 14 378-9 (1998)
# SOSUI: classification and secondary structure prediction system for membrane proteins.
#
# Mitaku S., Hirokawa T. Protein Eng. 11 (1999) Physicochemical factors for discriminating between soluble and membrane proteins: hydrophobicity of helical segments and protein length
#
# Mitaku S., Hirokawa T., and Tsuji T., Bioinformatics, 18 608-16 (2002)
# Amphiphilicity index of polar amino acids as an aid in the characterization of amino acid preference at membrane-water interfaces.



use strict;
use Getopt::Long;
use Data::Dumper;

my $input_file;
my $input_file_2;
GetOptions( "infile=s" => \$input_file, "infile2=s" => \$input_file_2);
if ((!defined $input_file) || (!defined $input_file_2)) {
	die "Usage $0 -infile SEQUENCES -infile2 PREDICTIONS\n";
}

my $output_file = "$input_file.sosui";
open( OUTFILE, ">$output_file") or
	die "Cannot open $output_file for writing\n";

open INFILE2, $input_file_2 or die $!;
my @input_lines_2 = <INFILE2>;

my %predictions;

my $reached_end = 0;
my $i = 0;
while ($reached_end == 0) {
	if ($i >= @input_lines_2) {
		$reached_end = 1;
	} else {
		my $input_line = $input_lines_2[$i];
		chomp($input_line);
		$input_line = trim($input_line);
		if ($input_line eq '') {
			$reached_end = 1;
		} else {
			my @bits = split(/ /, $input_line);
			my $this_fastaid = $bits[0];
			$this_fastaid = uc($this_fastaid);
			$i++;
			$input_line = $input_lines_2[$i];
			chomp($input_line);
			$input_line = trim($input_line);
			if ($input_line eq '') {
				$reached_end = 1;
			} else {
				if (substr($input_line,0,1) ne '>') {
					@bits = split(/REGION /, $input_line);
					my $this_tms = $bits[1] . ',';
					$predictions{$this_fastaid} = $this_tms;
					$i++;
				}
			}
		}
	}
}

open INFILE, $input_file or die $!;
my @input_lines = <INFILE>;

foreach my $input_line (@input_lines) {
	chomp($input_line);
	$input_line = trim($input_line);
	if ($input_line ne '') {
		my @bits = split(/:::::/, $input_line);
		my $seqid = $bits[0];
		my $fastaid = $bits[1];
		my $fastaid_key = uc($fastaid);
		my $aaseq = $bits[2];
		my $tms = '-';
		if (defined($predictions{$fastaid_key})) {
			$tms = $predictions{$fastaid_key};
			if ($tms eq '') {
				$tms = '-';
			}
		}
		my $output_line = $seqid . ':::::' . $fastaid . ':::::' . $aaseq . ':::::' . $tms . ':::::';
		print OUTFILE "$output_line\n";
	}
}
close OUTFILE;

sub trim {
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}




