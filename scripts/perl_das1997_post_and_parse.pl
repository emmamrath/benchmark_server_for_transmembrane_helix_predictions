#!/usr/bin/perl   -w

# perl perl_das1997_post_and_parse.pl -infile identity_30_graded.1line -cutoff 1.7 -rmvggg
# perl perl_das1997_post_and_parse.pl -infile scrap.fasta.1line -cutoff 2.2 -rmvggg

# If the -rmvggg flag is present, then remove leading and trailing GGG, HHH and AAA before submitting to prediction tool.

# will produce output file xxx.das1997 (eg. identity_30_graded.1line.das1997)

# Sequences were submitted to the DAS server http://www.sbc.su.se/~miklos/DAS/
# to predict transmembrane helices.

# Example of infile :
# 1uaz_A:::::>1uaz_A:::::TAAVGADLLGDGRPETLWLGIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIGLTEVQVGSEMLDIYYARYADWLFTTPLLLLDLALLAKVDRVSIGTLVGVDALMIVTGLVGALSHTPLARYTWWLFSTICMIVVLYFLATSLRAAAKERGPEVASTFNTLTALVLVLWTAYPILWIIGTEGAGVVGLGIETLLFMVLDVTAKVGFGFILLRSRAILGDTEAPEPSAGAEASAAD:::::
# 3ddl_A:::::>3ddl_A:::::MLQELPTLTPGQYSLVFNMFSFTVATMTASFVFFVLARNNVAPKYRISMMVSALVVFIAGYHYFRITSSWEAAYALQNGMYQPTGELFNDAYRYVDWLLTVPLLTVELVLVMGLPKNERGPLAAKLGFLAALMIVLGYPGEVSENAALFGTRGLWGFLSTIPFVWILYILFTQLGDTIQRQSSRVSTLLGNARLLLLATWGFYPIAYMIPMAFPEAFPSNTPGTIVALQVGYTIADVLAKAGYGVLIYNIAKAKSEEEGFNVSEMVEPATASA:::::

# Please cite: M. Cserzo, E. Wallin, I. Simon, G. von Heijne and A. Elofsson: 
# Prediction of transmembrane alpha-helices in procariotic membrane proteins: the Dense Alignment Surface method; 
# Prot. Eng. vol. 10, no. 6, 673-676, 1997

# Your query is:
# MVGLTTLFWLGAIGMLVGTLAFAWAGRDAGSGERRYYVTLVGISGIAAVAYVVMALGVGW
# GPPGVALLTPTVDVALIVYLDLVTKVGFGFIALDAAATLRAEHGESLAGVDTDAPAVAD
# 
# Potential transmembrane segments
# Start	Stop	Length	~	Cutoff
# 7	21	15	~	1.7
# 9	19	11	~	2.2
# 40	59	20	~	1.7
# 
# The DAS curve for your query:

# Done
# <hr><pre><em>Your query is:</em>
# MVGLTTLFWLGAIGMLVGTLAFAWAGRDAGSGERRYYVTLVGISGIAAVAYVVMALGVGW
# GPPGVALLTPTVDVALIVYLDLVTKVGFGFIALDAAATLRAEHGESLAGVDTDAPAVAD
# <hr>
# Potential transmembrane segments
# Start	Stop	Length	~	Cutoff
# 7	21	15	~	1.7
# 9	19	11	~	2.2
# 40	59	20	~	1.7
# 42	58	17	~	2.2
# 77	91	15	~	1.7
# 80	90	11	~	2.2
# 101	114	14	~	1.7
# 102	111	10	~	2.2
# 129	142	14	~	1.7
# 131	141	11	~	2.2
# 166	181	16	~	1.7
# 168	179	12	~	2.2
# 191	204	14	~	1.7
# 194	202	9	~	2.2
# 206	206	1	~	1.7
# 208	211	4	~	1.7
# <hr>
# The DAS curve for your query: </pre>
# 
# <img src="das1997_result_files/25497.gif">



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
my $url_for_submit = 'http://www.sbc.su.se/~miklos/DAS/tmdas.cgi';

# global variables

my $max_num_retries_of_query1 = 10;
my $count_input_seqs = 0;
my $referer = 'perl program';

my $input_file;
my $input_cutoff;
my $input_rmvggg;
GetOptions( "infile=s" => \$input_file, "cutoff=s" => \$input_cutoff, "rmvggg" => \$input_rmvggg );
if (!defined $input_file) {
	die "Usage $0 -infile FASTAfile -pat pattern\n";
}
if (!defined $input_cutoff) {
	$input_cutoff = 1.7;
}

# Open output files

my $output_file = "$input_file.das1997";
my $output_file_2 = "$input_file.das1997_html";
my $redo_file = "$input_file.das1997_redo";

open( OUTFILE, ">$output_file") or
	die "Cannot open $output_file for writing\n";
if ($create_outfile2 == 1) {
	open( OUTFILE2, ">$output_file_2") or
		die "Cannot open $output_file for writing\n";
}
open( REDOFILE, ">$redo_file") or
	die "Cannot open $output_file for writing\n";

open INFILE, $input_file or die $!;
my @input_lines = <INFILE>;

my $seqid;
my $fastaid;
my $aaseq;

foreach my $input_line (@input_lines) {
	chomp($input_line);
	$input_line = trim($input_line);
	if ($input_line ne '') {
		my @bits = split(/:::::/, $input_line);
		$seqid = $bits[0];
		$fastaid = $bits[1];
		$aaseq = $bits[2];
		process_1_fasta_seq();
	}
}

print OUTFILE "\n";
close OUTFILE;
if ($create_outfile2 == 1) {
	print OUTFILE2 "\n";
	close OUTFILE2;
}
close REDOFILE;

sub process_1_fasta_seq {

	my $input_aaseq_fix = $aaseq;
	$input_aaseq_fix =~ s/X/G/g; # convert invalid residue to valid, make sure program doesn't ignore residue and count residue positions incorrectly 

	my %seq_and_index_array;
	my $seq_and_index = \%seq_and_index_array;
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
	print         "$time_str processing input number $count_input_seqs, sequence $seqid...\n";
	# print OUTFILE "$time_str Processing $seqid\n";
	# print OUTFILE "$input_aaseq\n";
	if ($create_outfile2 == 1) {
		print OUTFILE2 "$time_str Processing $seqid\n\n";
	}

	my $ua = LWP::UserAgent->new;
	my $newagent = 'test search program (' . $ua->agent . ')';
	$ua->agent($newagent);

	my $request = POST ( $url_for_submit,
		Content => [ 	"QUERY .." => $seq_and_index->{'seq'}
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

		my @lines = split(/\n/, $content);

		my $output_info = '';

		my $in_tms_area = 0;
		my $finished_tms_area = 0;

		foreach my $line (@lines) {

			if ($finished_tms_area == 0) {
				if ($in_tms_area == 0) {
					if (substr($line, 0, 10) eq "Start\tStop") {
						$in_tms_area = 1;
					}
				} else { # $in_tms_area == 1
					if (substr($line, 0, 10) eq "\<hr\>") {
						$in_tms_area = 0;
						$finished_tms_area = 1;
					} else {
						my @bits = split(/\t/, $line);
						my $start = $bits[0] + $seq_and_index->{'index'};
						my $stop = $bits[1] + $seq_and_index->{'index'};
						my $cutoff = $bits[4];
						# loose cutoff = 1.7
						# strict cutoff = 2.2
						if ($cutoff == $input_cutoff) {
							$output_info .= $start . '-' . $stop . ',';
						}
					}
				}
			}
		}

		if ($output_info eq '') {
			$output_info = '-';
		}
		my $output_line = $seqid . ':::::' . $fastaid . ':::::' . $aaseq . ':::::' . $output_info . ':::::';

		print OUTFILE "$output_line\n";

	} else {
		print OUTFILE "Had $max_num_retries_of_query1 timeouts - got no results.\n";
		print REDOFILE "\>$seqid\n$aaseq\n\n";
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



sub trim {
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}




