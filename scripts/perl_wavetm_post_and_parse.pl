#!/usr/bin/perl   -w

# perl perl_wavetm_post_and_parse.pl -infile test1.fasta
# perl perl_wavetm_post_and_parse.pl -infile test1.fasta -rmvggg

# If the -rmvggg flag is present, then remove leading and trailing GGG, HHH and AAA before submitting to prediction tool.

# will produce output file xxx.wavetm (eg. scrap.fasta.wavetm)
# will produce output file xxx.wavetm_html (eg. scrap.fasta.wavetm_html)
# will produce output file xxx.wavetm_redo (eg. scrap.fasta.wavetm_redo)

# The input to this program is a file containing multiple fasta sequences.
# Each fasta sequence is submitted to the waveTM server at http://athina.biol.uoa.gr/bioinformatics/waveTM/input.html
# The REDO file looks just like the input file, and contains the sequences that always got rc 500 timeout,




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
my $url_for_submit = 'http://athina.biol.uoa.gr/cgi-bin/WAVELET/wave.pl';
#my $url_for_submit = 'http://127.0.0.1/benchmark_new/perlpgm3.pl';

# global variables

my $max_num_retries_of_query1 = 1;
my $count_input_seqs = 0;
my $output_info;
my $referer = 'http://athina.biol.uoa.gr/bioinformatics/waveTM/input.html';

my $fasta_file;
my $input_rmvggg;
GetOptions( "infile=s" => \$fasta_file, "rmvggg" => \$input_rmvggg );
if (!defined $fasta_file) {
	die "Usage $0 -infile FASTAfile\n";
}

# Open output files

my $output_file = "$fasta_file.wavetm";
my $output_file_link = "$fasta_file.wavetm_link";
my $output_file_2 = "$fasta_file.wavetm_html";
my $redo_file = "$fasta_file.wavetm_redo";

open( OUTFILE, ">$output_file") or
	die "Cannot open $output_file for writing\n";
open( OUTFILE_LINK, ">$output_file_link") or
	die "Cannot open $output_file_link for writing\n";
if ($create_outfile2 == 1) {
	open( OUTFILE2, ">$output_file_2") or
		die "Cannot open $output_file for writing\n";
}
open( REDOFILE, ">$redo_file") or
	die "Cannot open $output_file for writing\n";

my $in  = Bio::SeqIO->new(-file => $fasta_file , '-format' => 'Fasta');

while ( my $seq = $in->next_seq() ) {
	process_1_fasta_seq($seq);
}

print OUTFILE "\n";
close OUTFILE;
print OUTFILE_LINK "\n";
close OUTFILE_LINK;
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
	print         "$time_str processing input number $count_input_seqs, sequence $input_seqid...\n";
	# print OUTFILE "$time_str Processing $input_seqid\n";
	# print OUTFILE "$input_aaseq\n";
	if ($create_outfile2 == 1) {
		print OUTFILE2 "$time_str Processing $input_seqid\n\n";
	}

	my $ua = LWP::UserAgent->new;
	my $newagent = 'Mozilla/5.0 (X11; Ubuntu; Linux i686; rv:10.0.2) Gecko/20100101 Firefox/10.0.2';
	$ua->agent($newagent);

	my $input_fastaid_aaseq = "$input_fastaid\n$input_aaseq";

	#$seq_and_index->{'seq'} = 'MVEPLLSGIVLGLIVVTLAGLFYAAYKQYKRPNELGGGGGMV';

	my $seq_for_posting = $seq_and_index->{'seq'};

	my $seq_is_long_enough_for_wavetm = 0;
	while ($seq_is_long_enough_for_wavetm == 0) {
		if ( length($seq_for_posting) < 42 ) {
			$seq_for_posting .= $seq_for_posting;
		} else {
			$seq_is_long_enough_for_wavetm = 1;
		}
	}

	#my $seq_for_posting = $seq_and_index->{'seq'} . "\n";
	#my $seq_for_posting = "MWELRSASFWRAIFAEFFATLFYVFFGLGASLRWAPGPLHVLQVALAFGLALATLVQAVG\nHISGAHVNPAVTFAFLVGSQMSLLRAICYVVAQLLGAVAGAAVLYSVTPPAVRGNLALNT\nLHPGVSVGQATIVEIFLTLQFVLCIFATYDERRNGRLGSVALAVGFSLTLGHLFGMYYTG\nAGMNPARSFAPAILTRNFTNHWVYWVGPVIGAGLGSLLYDFLLFPRLKSVSERLSILKGT\nRPSESNGQPEVTGEPVELKTQAL";
#	my $seq_for_posting = "MWELRSASFWRAIFAEFFATLFYVFFGLGASLRWAPGPLHVLQVALAFGLALATLVQAVGHISGAHVNPAVTFAFLVGSQMSLLRAICYVVAQLLGAVAGAAVLYSVTPPAVRGNLALNTLHPGVSVGQATIVEIFLTLQFVLCIFATYDERRNGRLGSVALAVGFSLTLGHLFGMYYTGAGMNPARSFAPAILTRNFTNHWVYWVGPVIGAGLGSLLYDFLLFPRLKSVSERLSILKGT\nRPSESNGQPEVTGEPVELKTQAL";
	#my $seq_for_posting = "MWELRSASFWRAIFAEFFATLFYVFFGLGASLRWAPGPLHVLQVALAFGLALATLVQAVGHISGAHVNPAVTFAFLVGSQMSLLRAICYVVAQLLGAVAGAAVLYSVTPPAVRGNLALNT\nLHPGVSVGQATIVEIFLTLQFVLCIFATYDERRNGRLGSVALAVGFSLTLGHLFGMYYTGAGMNPARSFAPAILTRNFTNHWVYWVGPVIGAGLGSLLYDFLLFPRLKSVSERLSILKGTRPSESNGQPEVTGEPVELKTQAL";
	#my $seq_for_posting = "MWELRSASFWRAIFAEFFATLFYVFFGLGASLRWAPGPLHVLQVALAFGLALATLVQAVGHISGAHVNPAVTFAFLVGSQMSLLRAICYVVAQLLGAVAGAAVLYSVTPPAVRGNLALNT" . "\n" . "LHPGVSVGQATIVEIFLTLQFVLCIFATYDERRNGRLGSVALAVGFSLTLGHLFGMYYTGAGMNPARSFAPAILTRNFTNHWVYWVGPVIGAGLGSLLYDFLLFPRLKSVSERLSILKGTRPSESNGQPEVTGEPVELKTQAL";

	my $half_length = int( length($seq_for_posting) / 2 );
	$seq_for_posting = substr( $seq_for_posting, 0, $half_length ) . "\n" . substr( $seq_for_posting, $half_length );

	my $request = POST ( $url_for_submit,
		Content => [ 	seq_name => $input_seqid,
				#seq_name => '2zt9_G',
				seq => $seq_for_posting,
				#seq => 'MVEPLLSGIVLGLIVVTLAGLFYAAYKQYKRPNELGG',
				#seq => "MWELRSASFWRAIFAEFFATLFYVFFGLGASLRWAPGPLHVLQVALAFGLALATLVQAVG\nHISGAHVNPAVTFAFLVGSQMSLLRAICYVVAQLLGAVAGAAVLYSVTPPAVRGNLALNT\nLHPGVSVGQATIVEIFLTLQFVLCIFATYDERRNGRLGSVALAVGFSLTLGHLFGMYYTG\nAGMNPARSFAPAILTRNFTNHWVYWVGPVIGAGLGSLLYDFLLFPRLKSVSERLSILKGT\nRPSESNGQPEVTGEPVELKTQAL",
				show_obsTM => 'no_obs',
				obstm => '25-44, 89-110'
				#'' => 'Execute'
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

		# if ($create_outfile2 == 1) {
		#	print OUTFILE2 "\n...Printing out \$response\n\n";
		#	print OUTFILE2 Dumper($response);
		#	print OUTFILE2 "\n";
		# }

		# <?xml version="1.0" encoding="utf-8"?>
		# <!DOCTYPE html
		# 	PUBLIC "-//W3C//DTD XHTML Basic 1.0//EN"
		# 	"http://www.w3.org/TR/xhtml-basic/xhtml-basic10.dtd">
		# <html xmlns="http://www.w3.org/1999/xhtml" lang="en-US"><head><title>waveTM Results</title>
		# <style type="text/css">
		# <!--/* <![CDATA[ */
		# body {margin-left:10%;margin-right:10%;} h2 {margin-top:1em; margin-bottom:3em; margin-left:-4%}
		# 
		# /* ]]> */-->
		# </style>
		# </head><body bgcolor="#FFF0C0"><p><h2><img src='/WAVELET/waveTM.gif' 
		# align='left'>Prediction of 
		# transmembrane segments in proteins</h2></p><H4>RUNNING BASIC CHECKS...</H4><p /><H4>ok</H4><p /><H4>RUNNING WAVELET DENOISING OF THE HYDROPATHY SIGNAL</H4> (please wait ...)<p /><H4>ok</H4><p /><H4>RUNNING THE DYNAMIC PROGRAMMING ALGORITHM</H4> (please wait ...)<p /><H4>ok</H4><p />If your browser doesnt take you, go to  <a href=/waveTM/temp/65119.html>page</a><meta http-equiv='Refresh' content='0; URL=/bioinformatics/waveTM/temp/65119.html'>
		# 
		# if ($response->message eq '500') {
		# 	$got_timeout = 1;
		# } else {
		# 	$got_timeout = 0;
		# }

		my @lines = split(/\n/, $content);

		my $output_info = '';

		my $seen_section_3 = 0;
		my $in_tms_area = 0;
		my $finished_tms_area = 0;

		my $results_link = '-';

		foreach my $line (@lines) {
			#print "$line\n\n";

			my $find_result_link = index( $line, 'If your browser doesnt take you, go to' );
			if ($find_result_link > -1) {
				my @bits = split( /If your browser doesnt take you, go to  \<a href\=/, $line );
				my $bit1 = $bits[1];
				my @bits2 = split( /\>page\<\/a\>/, $bit1 );
				my $link_number = $bits2[0];
				$results_link = 'http://athina.biol.uoa.gr/bioinformatics' . $link_number;
			}
			#print OUTFILE_LINK "$line\n";
		}

		my $output_line_link = $input_seqid . ':::::' . $results_link . ':::::';

		print OUTFILE_LINK "$output_line_link\n";
		close OUTFILE_LINK;
		open( OUTFILE_LINK, ">>$output_file_link") or
			die "Cannot open $output_file_link for writing\n";

		if ($results_link ne '-') {

			my $request2 = HTTP::Request->new( GET => $results_link );
			my $response2 = $ua->request($request2);

			# <?xml version="1.0" encoding="utf-8"?>
			# <!DOCTYPE html
			# 	PUBLIC "-//W3C//DTD XHTML Basic 1.0//EN"
			# 	"http://www.w3.org/TR/xhtml-basic/xhtml-basic10.dtd">
			# <html xmlns="http://www.w3.org/1999/xhtml" lang="en-US"><head><title>waveTM Results</title>
			# <style type="text/css">
			# <!--/* <![CDATA[ */
			# body {margin-left:10%;margin-right:10%;} h2 {margin-top:1em; margin-bottom:3em; margin-left:-4%}
			# 
			# /* ]]> */-->
			# </style>
			# </head><body bgcolor="#FFF0C0"><p><h2><img src='/WAVELET/waveTM.gif' 
			# align='left'>Prediction of transmembrane segments in proteins</h2></p><H4>PREDICTED TRANSMEMBRANE SEGMENTS</H4><table border=0 cellpadding=0 width=70% align=center><tr><td><table border=1 cellpadding=5><tr><th bgcolor='#FFCC99' align=center># TM</th><th bgcolor='#FFCC99' align=center>Begin</th><th bgcolor='#FFCC99' align=center>End</th></tr><tr><td bgcolor='#FFEBC1' align=center>1</td><td bgcolor='#FFEBC1' align=center>4</td><td bgcolor='#FFEBC1' align=center>26</td></tr><tr><td bgcolor='#FFEBC1' align=center>2</td><td bgcolor='#FFEBC1' align=center>41</td><td bgcolor='#FFEBC1' align=center>63</td></tr></table></td></tr></table>
			#  <H4>SEQUENCE: 2zt9_G</H4><pre><b>              1         2         3         4         5         6
			#      123456789012345678901234567890123456789012345678901234567890
			# 	 .    |    .    |    .    |    .    |    .    |    .    |
			# <u>0000 </u>MVE<font color="blue">PLLSGIVLGLIVVTLAGLFYAAY</font>KQYKRPNELGGMVE<font color="blue">PLLSGIVLGLIVVTLAGLFY
			# <font color="black"><u>0060 </u></font>AAY</font>KQYKRPNELGG</u>
			# </b></pre><H4><I>Predicted transmembrane regions are indicated in <font color="blue">blue</font><p>The plot is <a href=/bioinformatics/waveTM/temp/65140.ps>here</a></p><H3>Topology prediction for your sequence is provided, executing: <BR><BR> <a href="http://athina.biol.uoa.gr/cgi-bin/orienTM/orienTM?seq=ID+++Prediction+by+WaveTM%3B%0D%0AAC+++Prediction+by+WaveTM%3B%0D%0AFT+++TRANSMEM++++++4+++++26+++++++PREDICTION%3B%0D%0AFT+++TRANSMEM+++++41+++++63+++++++PREDICTION%3B%0D%0ASQ%0D%0A+++++MVEPLLSGIVLGLIVVTLAGLFYAAYKQYKRPNELGGMVEPLLSGIVLGLIVVTLAGLFY%0D%0A+++++AAYKQYKRPNELGG%0D%0A%2F%2F%0D%0A"><img src='http://athina.biol.uoa.gr/orienTM/orienTM.gif'></a><p>Back to <a href="http://biophysics.biol.uoa.gr">Biophysics Lab</a></p>

			if ($response2->is_success) {
				my $http_result2 = $response2->content;
				my @lines2 = split(/\n/, $http_result2);

				my $tmh_string = '';

				foreach my $line2 (@lines2) {
					#print "$line2\n";

					my $find_tmh_heading = index( $line2, 'Prediction of transmembrane segments in proteins' );
					if ($find_tmh_heading > -1) {
						my @bits = split( /End\<\/th\>\<\/tr\>/, $line2 );
						my $string_with_tmh_rows = $bits[1];
						my @rows = split( /\<tr\>/, $string_with_tmh_rows );
						foreach my $this_row (@rows) {
							$this_row = trim($this_row);
							if ($this_row ne '') {
								#print "$this_row\n";
								my @cells = split( /\<td bgcolor=\'\#FFEBC1\' align=center\>/, $this_row );
								my $this_from_plus_extra = $cells[2];
								my $this_to_plus_extra = $cells[3];
								my @bits2 = split( /\<\/td\>/, $this_from_plus_extra );
								my $this_from = $bits2[0];
								my @bits3 = split( /\<\/td\>/, $this_to_plus_extra );
								my $this_to = $bits3[0];
								$tmh_string .= $this_from . '-' . $this_to . ',';
								#print "$tmh_string\n\n";
							}
						}
					}
				}

				my $output_line = $input_seqid . ':::::' . $tmh_string . ':::::';

				print OUTFILE "$output_line\n";
				close OUTFILE;
				open( OUTFILE, ">>$output_file") or
					die "Cannot open $output_file for writing\n";
			}
		}

	} else {
		print OUTFILE "$input_seqid had $max_num_retries_of_query1 timeouts - got no results.\n";
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



sub trim {
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}

