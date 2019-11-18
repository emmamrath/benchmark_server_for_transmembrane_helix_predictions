#!/usr/bin/perl   -w

# perl perl_hmmtm_post_and_save.pl -infile test1.fasta

# If the -rmvggg flag is present, then remove leading and trailing GGG, HHH and AAA before submitting to prediction tool.

# will produce output file xxx.hmmtm_output (eg. scrap.fasta.hmmtm_output)
# will produce output file xxx.hmmtm_html (eg. scrap.fasta.hmmtm_html)
# will produce output file xxx.hmmtm_redo (eg. scrap.fasta.hmmtm_redo)

# The input to this program is a file containing multiple fasta sequences.
# Each fasta sequence is submitted to the HMM-TM server at http://bioinformatics.biol.uoa.gr/HMM-TM/input.jsp
# The REDO file looks just like the input file, and contains the sequences that always got rc 500 timeout,

# If you are going to use these results in your work, please cite:
# G.E Tusnády and I. Simon (1998) 
# J. Mol. Biol.283, 489-506.
# G.E Tusnády and I. Simon (2001) 
# Bioinformatics17, 849-850.



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
my $url_for_submit = 'http://bioinformatics.biol.uoa.gr/HMM-TM/process.jsp';

# global variables

my $max_num_retries_of_query1 = 10;
my $count_input_seqs = 0;
my $output_info;
my $referer = 'Client User Agent';

my $fasta_file;
my $input_rmvggg;
GetOptions( "infile=s" => \$fasta_file, "rmvggg" => \$input_rmvggg );
if (!defined $fasta_file) {
	die "Usage $0 -infile FASTAfile -pat pattern\n";
}

# Open output files

my $output_file = "$fasta_file.hmmtm_output";
my $output_file_2 = "$fasta_file.hmmtm_html";
my $redo_file = "$fasta_file.hmmtm_redo";

open( OUTFILE, ">$output_file") or
	die "Cannot open $output_file for writing\n";
if ($create_outfile2 == 1) {
	open( OUTFILE2, ">$output_file_2") or
		die "Cannot open $output_file for writing\n";
}
open( REDOFILE, ">$redo_file") or
	die "Cannot open $output_file for writing\n";

my $in  = Bio::SeqIO->new(-file => $fasta_file , '-format' => 'Fasta');

my %seq_and_index_array;
my $seq_and_index = \%seq_and_index_array;

while ( my $seq = $in->next_seq() ) {
	process_1_fasta_seq($seq);
}

print OUTFILE "\n";
close OUTFILE;
if ($create_outfile2 == 1) {
	print OUTFILE2 "\n";
	close OUTFILE2;
}
close REDOFILE;

sub process_1_fasta_seq {

	my $seq = shift;
	# print "process_1_fasta_seq: Sequence ",$seq->id," first 10 bases ",$seq->subseq(1,10),"\n";
	# print Dumper($seq);
	my $input_aaseq = $seq->seq;
	my $input_seqid = $seq->id . '';
	my $input_fastaid = '>' . $input_seqid;
	my $input_seq_desc = '';
	if (defined($seq->desc)) {
		$input_seq_desc = $seq->desc . '';
	}

	my $input_aaseq_fix = $input_aaseq;
	$input_aaseq_fix =~ s/X/G/g; # convert invalid residue to valid, make sure program doesn't ignore residue and count residue positions incorrectly 

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
	print         "$time_str processing input number $count_input_seqs, sequence $input_seqid...\n";
	# print OUTFILE "$time_str Processing $input_seqid\n";
	# print OUTFILE "$input_aaseq\n";
	if ($create_outfile2 == 1) {
		print OUTFILE2 "$time_str Processing $input_seqid\n\n";
	}

	my $ua = LWP::UserAgent->new;
	#my $newagent = 'test search program (' . $ua->agent . ')';
	my $newagent = 'Client User Agent';
	$ua->agent($newagent);

	my $input_fastaid_aaseq = "$input_fastaid\n$input_aaseq";

	my $request = POST ( $url_for_submit,
		Content => [ 	'wait.target' => 'process.jsp',
				'formInput' => $input_fastaid_aaseq,
				#'formInput' => $seq_and_index->{'seq'},
				#'method' => 0 # Viterbi
				#'method' => 1 # N-Best
				#'method' => 2 # Posterior
				'method' => 3, # Posterior Viterbi
				'formInput2' => '',
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

		$output_info = '';

		my $output_line = $input_seqid . ':::::';
		print OUTFILE "$output_line\n";

		# <table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>1&nbsp;&nbsp;&nbsp;&nbsp;4</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>5&nbsp;&nbsp;&nbsp;&nbsp;25</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>26&nbsp;&nbsp;&nbsp;37</span></td></tr></table></td>

		# <TABLE align="center">
		# <td class=col5><span class="seqcla2">
		# <table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>1&nbsp;&nbsp;&nbsp;&nbsp;9</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>10&nbsp;&nbsp;&nbsp;27</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>28&nbsp;&nbsp;&nbsp;46</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>47&nbsp;&nbsp;&nbsp;64</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>65&nbsp;&nbsp;&nbsp;74</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>75&nbsp;&nbsp;&nbsp;95</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>96&nbsp;&nbsp;&nbsp;102</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>103&nbsp;&nbsp;124</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>125&nbsp;&nbsp;144</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>145&nbsp;&nbsp;162</span></td></tr></table></td><td class=col5><span class="seqcla2">
		# <table><tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>163&nbsp;&nbsp;167</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>168&nbsp;&nbsp;186</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>187&nbsp;&nbsp;221</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>222&nbsp;&nbsp;239</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>240&nbsp;&nbsp;259</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>260&nbsp;&nbsp;278</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>279&nbsp;&nbsp;290</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>291&nbsp;&nbsp;309</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>310&nbsp;&nbsp;314</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>315&nbsp;&nbsp;331</span></td></tr></table></td><td class=col5><span class="seqcla2">
		# <table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>332&nbsp;&nbsp;346</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>347&nbsp;&nbsp;365</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>366&nbsp;&nbsp;380</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>381&nbsp;&nbsp;402</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>403&nbsp;&nbsp;417</span></td></tr></table></td>
		# </tr>
		# </table>

		foreach my $line (@lines) {

			#parse_1_html_search_output_line_of_1_fasta_seq($line);

			if ( index( $line, 'class="seqcla2"' ) > -1 ) {

				print OUTFILE "$line\n";
			}
		}

		#if ($output_info eq '') {
		#	$output_info = '-';
		#}

		#my $output_line = $input_seqid . ':::::' . $input_fastaid . ':::::' . $input_aaseq . ':::::' . $output_info . ':::::';

		#print OUTFILE "$output_line\n";

		close OUTFILE;
		open( OUTFILE, ">>$output_file") or
			die "Cannot open $output_file for writing\n";

	} else {
		print OUTFILE "Had $max_num_retries_of_query1 timeouts - got no results.\n";
		print REDOFILE "\>$input_seqid\n$input_aaseq\n\n";
	}
}

sub parse_1_html_search_output_line_of_1_fasta_seq {

	my $line = shift;

	if (substr($line, 0, 23) eq 'Transmembrane helices: ') {
		my $this_info = substr($line, 23);
		$this_info =~ s/ /\,/g;
		my @bits = split(/\,/, $this_info);
		foreach my $from_to (@bits) {
			my @bits2 = split(/\-/, $from_to);
			my $from = $bits2[0] + $seq_and_index->{'index'};
			my $to = $bits2[1] + $seq_and_index->{'index'};
			my $new_from_to = $from . '-' . $to . ',';
			$output_info .= $new_from_to;
		}
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

