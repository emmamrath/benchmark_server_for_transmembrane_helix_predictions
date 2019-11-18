#!/usr/bin/perl   -w

# perl perl_PHDThtm_PBIL_post_and_save.pl -infile 3bz1_D.fasta

# If the -rmvggg flag is present, then remove leading and trailing GGG, HHH and AAA before submitting to prediction tool.

# will produce output file xxx.PHDhtm_PBIL_output (eg. scrap.fasta.PHDThtm_PBIL_output)
# will produce output file xxx.PHDhtm_PBIL_html (eg. scrap.fasta.PHDThtm_PBIL_html)
# will produce output file xxx.PHDhtm_PBIL_redo (eg. scrap.fasta.PHDThtm_PBIL_redo)

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

my $output_file = "$fasta_file.PHDThtm_PBIL_output";
my $output_file_2 = "$fasta_file.PHDThtm_PBIL_html";
my $redo_file = "$fasta_file.PHDThtm_PBIL_redo";

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

		my $phdhtm_link = '';

		foreach my $line (@lines) {

			# line=Prediction result file (text): [<A HREF=/tmp/70952b26452d.phdhtm>PHD</A>]

			if ( index( $line, 'Prediction result file (text):' ) > -1 ) {
				my @bits = split(/\<A HREF\=/, $line);
				my $bit1 = $bits[1];
				my @bits2 = split(/\>PHD\<\/A\>/, $bit1);
				$phdhtm_link = 'http://npsa-pbil.ibcp.fr' . $bits2[0];
			}
		}

		#--- other        : predictions derived based on PHDhtm
		#--- PHDFhtm      : filtered prediction, i.e., too long HTM's are
		#---                split, too short ones are deleted
		#--- PHDRhtm      : refinement of neural network output 
		#--- PHDThtm      : topology prediction based on refined model
		#---                symbols used:
		#---             i: intra-cytoplasmic
		#---             T: transmembrane region
		#---             o: extra-cytoplasmic
		#--- 
		#--- PhdTopology REFINEMENT AND TOPOLOGY PREDICTION
		#		  ....,....1....,....2....,....3....,....4....,....5....,....6
		#	  AA      |MTIAIGRAPAERGWFDILDDWLKRDRFVFVGWSGILLFPCAYLALGGWLTGTTFVTSWYT|
		#	  PHD htm |                            HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH|
		#	  Rel htm |999999999999999999999999874024677888777777777777778877776531|
		# detail: 
		#	  prH htm |000000000000000000000000012467888999888888888888889988888765|
		#	  prL htm |999999999999999999999999987532111000111111111111110011111234|
		# subset: SUB htm |LLLLLLLLLLLLLLLLLLLLLLLLLL.....HHHHHHHHHHHHHHHHHHHHHHHHH....|
		#	  PHDRhtm |                              HHHHHHHHHHHHHHHHHHHHHHHHH     |
		#	  PHDThtm |iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiTTTTTTTTTTTTTTTTTTTTTTTTTooooo|
		#		  ....,....7....,....8....,....9....,....10...,....11...,....12
		#	  AA      |HGLASSYLEGCNFLTVAVSTPANSMGHSLLLLWGPEAQGDFTRWCQLGGLWTFIALHGAF|
		#	  PHD htm |             HHHHHHHHHHHHHH                     HHHHHHHHHHHH|
		#	  Rel htm |123344444431023566667766530245678999999999998763045677888877|
		# detail: 
		#	  prH htm |433322222234466788888888765322110000000000000113577888999988|
		#	  prL htm |566677777765533211111111234677889999999999999886422111000011|
		# subset: SUB htm |....................HH.........LLLLLLLLLLLLLLL......HHHHHHHH|
		#	  PHDRhtm |         HHHHHHHHHHHHHHHHHH                     HHHHHHHHHHHH|
		#	  PHDThtm |oooooooooTTTTTTTTTTTTTTTTTTiiiiiiiiiiiiiiiiiiiiiTTTTTTTTTTTT|
		#		  ....,....13...,....14...,....15...,....16...,....17...,....18
		#	  AA      |GLIGFMLRQFEIARLVGVRPYNAIAFSAPIAVFVSVFLIYPLGQSSWFFAPSFGVAAIFR|
		#	  PHD htm |HHHHHHH            HHHHHHHHHHHHHHHHHHHHHHHHHHHH   HHHHHHHHHH|
		#	  Rel htm |777764113566677542125677777777777777777776421000102567777877|
		# detail: 
		#	  prH htm |888887543211111223467888888888888888888888765554446788888988|
		# ...
		#		  ....,....31...,....32...,....33...,....34...,....35...,....36
		#	  AA      |QEIRAAEDPEFETFYTKNLLLNEGIRAWMAPQDQPHENFVFPEEVLPRGNAL|
		#	  PHD htm |                                                    |
		#	  Rel htm |9999999999999999999999999999999999999999999999999999|
		# detail: 
		#	  prH htm |0000000000000000000000000000000000000000000000000000|
		#	  prL htm |9999999999999999999999999999999999999999999999999999|
		# subset: SUB htm |LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL|
		#	  PHDRhtm |                                                    |
		#	  PHDThtm |oooooooooooooooooooooooooooooooooooooooooooooooooooo|
		#--- 
		#--- PhdTopology REFINEMENT AND TOPOLOGY PREDICTION END
		#--- 

		#print "phdhtm_link=$phdhtm_link\n";#####

		if ($phdhtm_link ne '') {

			my $request2 = HTTP::Request->new( GET => $phdhtm_link );
			my $response2 = $ua->request($request2);

			my $output_topology = '';

			if ($response2->is_success) {
				my $http_result2 = $response2->content;
				my @lines2 = split(/\n/, $http_result2);

				foreach my $line2 (@lines2) {

					if ( index( $line2, 'PHDThtm |' ) > -1 ) {
						my @bits3 = split(/\|/, $line2);
						my $bit3 = $bits3[1];
						$output_topology .= $bit3;
					}
				}
			}

			if ($output_topology eq '') {
				$output_topology = '-';
			}
			my $output_line = $input_seqid . ':::::' . $output_topology . ':::::';

			print OUTFILE "$output_line\n";

			close OUTFILE;
			open( OUTFILE, ">>$output_file") or
				die "Cannot open $output_file for rewriting\n";
		}

	} else {
		print OUTFILE "Had $max_num_retries_of_query1 timeouts - got no results.\n";
		print REDOFILE "\>$input_seqid\n$input_aaseq\n\n";
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

