#!/usr/bin/perl   -w

# perl perl_das2002_post_and_parse.pl -infile identity_30_graded.1line -rmvggg
# perl perl_das2002_post_and_parse.pl -infile scrap.fasta.1line

# If the -rmvggg flag is present, then remove leading and trailing GGG, HHH and AAA before submitting to prediction tool.

# will produce output file xxx.das2002 (eg. identity_30_graded.1line.das2002)

# Sequences were submitted to the DAS server http://mendel.imp.ac.at/DAS/
# (was http://mendel.imp.ac.at/sat/DAS/DAS.html )
# to predict transmembrane helices.

# Example of infile :
# 1uaz_A:::::>1uaz_A:::::TAAVGADLLGDGRPETLWLGIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIGLTEVQVGSEMLDIYYARYADWLFTTPLLLLDLALLAKVDRVSIGTLVGVDALMIVTGLVGALSHTPLARYTWWLFSTICMIVVLYFLATSLRAAAKERGPEVASTFNTLTALVLVLWTAYPILWIIGTEGAGVVGLGIETLLFMVLDVTAKVGFGFILLRSRAILGDTEAPEPSAGAEASAAD:::::
# 3ddl_A:::::>3ddl_A:::::MLQELPTLTPGQYSLVFNMFSFTVATMTASFVFFVLARNNVAPKYRISMMVSALVVFIAGYHYFRITSSWEAAYALQNGMYQPTGELFNDAYRYVDWLLTVPLLTVELVLVMGLPKNERGPLAAKLGFLAALMIVLGYPGEVSENAALFGTRGLWGFLSTIPFVWILYILFTQLGDTIQRQSSRVSTLLGNARLLLLATWGFYPIAYMIPMAFPEAFPSNTPGTIVALQVGYTIADVLAKAGYGVLIYNIAKAKSEEEGFNVSEMVEPATASA:::::

# Please cite in your relevant publications:
# "On filtering false positive transmembrane protein predictions": 
# Miklos Cserzo, Frank Eisenhaber, Birgit Eisenhaber, and Istvan Simon; 
# Protein Eng. 2002 15: 745-752.

# === Result of the prediction ===
# 
# &gt;1uaz_A
# # TMH:  7 Q: trusted
# @   27   4.887 core:   19 ..   34 1.097e-05
# @   63   3.104 core:   54 ..   69 5.936e-03
# @  100   3.914 core:   94 ..  106 3.406e-04
# @  124   3.878 core:  117 ..  131 3.859e-04
# @  153   5.657 core:  145 ..  160 7.247e-07
# @  186   4.936 core:  179 ..  197 9.229e-06
# @  213   3.597 core:  205 ..  224 1.042e-03
# 
# &lt;-------- end of list --------&gt;



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
my $url_for_submit = 'http://mendel.imp.ac.at/DAS/cgi-bin/das.cgi';
# was my $url_for_submit = 'http://mendel.imp.ac.at/sat/DAS/cgi-bin/das.cgi';

# global variables

my $max_num_retries_of_query1 = 10;
my $count_input_seqs = 0;
my $referer = 'perl program';

my $input_file;
my $input_rmvggg;
GetOptions( "infile=s" => \$input_file, "rmvggg" => \$input_rmvggg );
if (!defined $input_file) {
	die "Usage $0 -infile INFILE\n";
}

# Open output files

my $output_file = "$input_file.das2002";
my $output_file_2 = "$input_file.das2002_html";
my $redo_file = "$input_file.das2002_redo";

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

	# <b>Output format:</b> long
	# <input name="LEN" value="-l" type="radio">
	# short
	# <input name="LEN" value="-s" checked="checked" type="radio">
	# 
	# </p><p>
	# <b>Evaluation:</b> unconditional
	# <input name="CON" value="-u" type="radio">
	# trusted
	# <input name="CON" value="-t" checked="checked" type="radio">
	# 
	# </p><p>
	# <b>TM-library size:</b>
	# 8
	# <input name="LIB" value="08" checked="checked" type="radio">
	# 16
	# <input name="LIB" value="16" type="radio">
	# 24
	# <input name="LIB" value="24" type="radio">
	# 32
	# <input name="LIB" value="32" type="radio">

	my $input_field = "$fastaid\n$aaseq";
	my $request = POST ( $url_for_submit,
		Content => [ 	LEN => "-s",
				CON => "-t",
				LIB => "08",
				SEQ => $seq_and_index->{'seq'}
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

		foreach my $line (@lines) {

			if (substr($line, 0, 1) eq '@') {
				# @   27   4.887 core:   19 ..   34 1.097e-05
				# @   63   3.104 core:   54 ..   69 5.936e-03
				my @bits = split(/[\s]+/, $line);
				my $start = $bits[4] + $seq_and_index->{'index'};
				my $stop = $bits[6] + $seq_and_index->{'index'};
				$output_info .= $start . '-' . $stop . ',';
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




