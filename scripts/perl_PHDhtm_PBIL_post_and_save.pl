#!/usr/bin/perl   -w

# perl perl_PHDhtm_PBIL_post_and_save.pl -infile test1.fasta

# If the -rmvggg flag is present, then remove leading and trailing GGG, HHH and AAA before submitting to prediction tool.

# will produce output file xxx.PHDhtm_PBIL_output (eg. scrap.fasta.PHDhtm_PBIL_output)
# will produce output file xxx.PHDhtm_PBIL_html (eg. scrap.fasta.PHDhtm_PBIL_html)
# will produce output file xxx.PHDhtm_PBIL_redo (eg. scrap.fasta.PHDhtm_PBIL_redo)

# The input to this program is a file containing multiple fasta sequences.
# Each fasta sequence is submitted to the HMM-TM server at http://npsa-pbil.ibcp.fr/cgi-bin/npsa_automat.pl?page=/NPSA/npsa_htm.html
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
my $url_for_submit = 'http://npsa-pbil.ibcp.fr/cgi-bin/primanal_htm.pl';

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

my $output_file = "$fasta_file.PHDhtm_PBIL_output";
my $output_file_2 = "$fasta_file.PHDhtm_PBIL_html";
my $redo_file = "$fasta_file.PHDhtm_PBIL_redo";

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
		Content => [ 	'title' => $input_seqid,
				'notice' => $input_aaseq,
				'ali_width' => 9999
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

		#<table border="0" cellpadding="0"><tbody><tr><td valign="middle" width="47" align="center"><img src="2zt9_G_PHDhtm_results_files/mechanism.gif" width="47" height="35" align="middle"></td><td valign="middle" width="600" align="center"><font size="-1"> Job <b>PHDHTM</b> (ID: 3c1c4a3b0498) is running on <font color="red"><b>NPS@</b></font> server (started on 20120323-075008).<br>Results will be shown below. <b>Please wait and don't go back.</b></font></td></tr></tbody></table>
		#<hr>
		#<font size="-1"><b>In your publication cite : </b><br>NPS@: Network Protein Sequence Analysis<br>TIBS 2000 March Vol. 25, No 3 [291]:147-150<br>Combet C., Blanchet C., Geourjon C. and Deléage G.<br></font><hr><h2>PHD transmembrane #helix prediction result for : 2zt9_G </h2>
		#<a href="http://npsa-pbil.ibcp.fr/NPSA/npsa_references.html#phdhtm">Abstract</a> <font size="-1">Rost B, Casadio R, Fariselli P, Sander C : Transmembrane helices predicted at 95% accuracy. 
		#					Protein Sci. 1995 Mar;4(3):521-33.</font><pre><code>
		#        10        20        30
		#         |         |         |
		#MVEPLLSGIVLGLIVVTLAGLFYAAYKQYKRPNELGG
		#    <font color="red">HHHHHHHHHHHHHHHHHHHH</font>             
		#
		#
		#
		#</code></pre>
		#<p>
		#Prediction result file (text): [<a href="http://npsa-pbil.ibcp.fr/tmp/3c1c4a3b0498.phdhtm">PHD</a>]
		#<br>
		#Intermediate result file (text): [<a href="http://npsa-pbil.ibcp.fr/tmp/3c1c4a3b0498.blastpsecpred">BLASTP on NRPROT</a>] [<a href="http://npsa-pbil.ibcp.fr/tmp/3c1c4a3b0498.cluoutphdpred">CLUSTALW</a>]
		#</p>
		#<hr><i><b>User :</b> public@210.84.16.171. <b>Last modification time :</b> Fri Mar 23 07:50:12 2012. <b>Current time :</b> Fri Mar 23 07:50:12 2012</i>

		my $in_output = 0;
		my $output_finished = 0;

		foreach my $line (@lines) {

			my $line_uc = uc($line);

			#parse_1_html_search_output_line_of_1_fasta_seq($line);

			if ($output_finished == 0) {
				if ( $in_output == 1) {
					if ( index( $line_uc, '<HR><I><B>USER' ) > -1 ) {
						$output_finished = 1;
						$in_output = 0;
					}
				} elsif ( index( $line_uc, '<PRE><CODE>' ) > -1 ) {
					$in_output = 1;
				}
			}
			#print "output_finished=$output_finished, in_output=$in_output\n$line\n\n";

			if ($in_output == 1) {
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

