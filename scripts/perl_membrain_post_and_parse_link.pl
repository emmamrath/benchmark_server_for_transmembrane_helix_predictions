#!/usr/bin/perl   -w

# perl perl_membrain_post_and_parse_link.pl -infile test1.fasta
# perl perl_membrain_post_and_parse_link.pl -infile test1.fasta -rmvggg

# If the -rmvggg flag is present, then remove leading and trailing GGG, HHH and AAA before submitting to prediction tool.

# will produce output file xxx.membrain (eg. scrap.fasta.membrain)
# will produce output file xxx.membrain_html (eg. scrap.fasta.membrain_html)
# will produce output file xxx.membrain_redo (eg. scrap.fasta.membrain_redo)

# The input to this program is a file containing multiple fasta sequences.
# Each fasta sequence is submitted to the MemBrain server at http://www.csbio.sjtu.edu.cn/bioinf/MemBrain/
# The REDO file looks just like the input file, and contains the sequences that always got rc 500 timeout,

# MemBrain demonstrates an overall improvement of about 20% in prediction accuracy,
# particularly, in predicting the ends of TMHs and TMHs that are shorter than 15 residues. 
# It also has the capability to detect N-terminal signal peptides. 

# Hongbin Shen, James J. Chou
# MemBrain: Improving the Accuracy of Predicting Transmembrane Helices
# PLoS ONE, 2008, 6: e2399. 



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
my $url_for_submit = 'http://www.csbio.sjtu.edu.cn/cgi-bin/MEMBRAIN.cgi';

# global variables

my $max_num_retries_of_query1 = 10;
my $count_input_seqs = 0;
my $output_info;
my $referer = 'Client User Agent';

my $fasta_file;
my $input_rmvggg;
GetOptions( "infile=s" => \$fasta_file, "rmvggg" => \$input_rmvggg );
if (!defined $fasta_file) {
	die "Usage $0 -infile FASTAfile\n";
}

# Open output files

my $output_file = "$fasta_file.membrain";
my $output_file_link = "$fasta_file.membrain_link";
my $output_file_2 = "$fasta_file.membrain_html";
my $redo_file = "$fasta_file.membrain_redo";

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
	my $newagent = 'Client User Agent (' . $ua->agent . ')';
	$ua->agent($newagent);

	my $input_fastaid_aaseq = "$input_fastaid\n$input_aaseq";

	my $request = POST ( $url_for_submit,
		Content => [ 	mode => 'string',
				S1 => $seq_and_index->{'seq'},
				R2 => 'SignalNO', # Checked : I know there is NO N-terminal signal peptide
				email => 'patricia@kiddiesgames.com'
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

		# <html><head>
		# <meta http-equiv="content-type" content="text/html; charset=GB2312">
		# <title>MemBrain: Improving the accuracy of predicting transmembrane helices</title><link rel="stylesheet" type="text/css" href="2b6o_A_membrain_link_files/css.html"></head>
		#
		# <body>
		# <div align="center"><table style="border-right: 1px solid rgb(128, 128, 128); border-left: 1px solid rgb(128, 128, 128);" id="table1" width="730" bgcolor="#f0f8ff" cellpadding="0" cellspacing="0"><tbody><tr><td style="border-bottom: 1px solid rgb(161, 161, 161);" bgcolor="#dff0ff" height="80"><table id="table2"><tbody><tr><td><table width="730"><tbody><tr><td><font size="5pt" face="times new roman"><b>MemBrain</b>:</font><font size="5pt" face="times new roman"> Improving the accuracy of predicting transmembrane helices</font></td></tr></tbody></table></td></tr><tr><td valign="bottom" height="40"><p align="left"> &nbsp;<font size="4pt"> <font face="times new roman">|</font> <a target="_blank" href="http://www.csbio.sjtu.edu.cn/bioinf/MemBrain/Readme.htm"><font color="blue" face="times new roman"><u>Read Me</u></font></a> &nbsp;<font face="times new roman">|</font> &nbsp;<a target="_blank" href="http://www.csbio.sjtu.edu.cn/bioinf/MemBrain/Data.htm"><font color="blue" face="times new roman"><u>Data</u></font></a> &nbsp;<font face="times new roman">|</font> &nbsp;<a target="_blank" href="http://www.csbio.sjtu.edu.cn/bioinf/MemBrain/Citation.htm"><font color="blue" face="times new roman"><u>Citation</u></font></a> &nbsp;</font>&nbsp; </p></td></tr></tbody></table></td></tr><tr><td>&nbsp;</td></tr><tr><td><font size="+1" face="times new roman">Thank you for using MemBrain and we will send the results to your email.<br> You can also click <font color="blue"><a href="http://www.csbio.sjtu.edu.cn/temp/11158/Result.txt" target="_new">http://www.csbio.sjtu.edu.cn/temp/11158/Result.txt</a></font> for the prediction results and click <font color="blue"><a href="http://www.csbio.sjtu.edu.cn/temp/11158/propensity.gif" target="_new">http://www.csbio.sjtu.edu.cn/temp/11158/propensity.gif</a></font> to see the plot of the TMH propensity of your query sequence when the job is finished. <br>Thank you for your patience!</font></td></tr><tr><td>&nbsp;</td></tr><tr><td align="center"><a href="http://www.csbio.sjtu.edu.cn/bioinf/MemBrain/"><font size="5pt" color="blue" face="times new roman">Home Page</font></a></td></tr><tr><td>&nbsp;</td></tr><tr><td style="border-top: 1px solid rgb(128, 128, 128);" align="center" bgcolor="#dff0ff" height="50"><font size="4pt" face="times new roman">Contact @ <a href="mailto:hbshen@sjtu.edu.cn">Hong-Bin</a></font></td></tr></tbody></table></div>
		# 
		# </body></html>

		my $content = $response->content;

		# if ($create_outfile2 == 1) {
		#	print OUTFILE2 "\n...Printing out \$content\n\n";
		#	print OUTFILE2 "$content\n\n";
		# }

		my @lines = split(/\n/, $content);

		my $output_info = '';

		my $seen_section_3 = 0;
		my $in_tms_area = 0;
		my $finished_tms_area = 0;

		my $results_link = '-';

		foreach my $line (@lines) {

			my $find_result_link = index( $line, 'Result.txt' );
			if ($find_result_link > -1) {
				my @bits = split( /http\:\/\/www.csbio.sjtu.edu.cn\/temp\//, $line );
				my $bit1 = $bits[1];
				my @bits2 = split( /\/Result\.txt/, $bit1 );
				my $link_number = $bits2[0];
				$results_link = 'http://www.csbio.sjtu.edu.cn/temp/' . $link_number . '/Result.txt';
			}
		}

		my $output_line_link = $input_seqid . ':::::' . $results_link . ':::::';

		print OUTFILE_LINK "$output_line_link\n";
		close OUTFILE_LINK;
		open( OUTFILE_LINK, ">>$output_file_link") or
			die "Cannot open $output_file_link for writing\n";

		if ($results_link ne '-') {

			my $request2 = HTTP::Request->new( GET => $results_link );
			my $response2 = $ua->request($request2);

			# ***************************************************************
			# 		MemBrain results notification                 
			# ***************************************************************
			# 
			# Thank you for using MemBrain. Predicted results are shown below and for any
			# question, please contact with Dr.James Chou at: james_chou@hms.harvard.edu or
			# Dr.Hong-Bin Shen at: hbshen@sjtu.edu.cn
			# Group webpage: http://chou.med.harvard.edu/ or http://www.csbio.sjtu.edu.cn
			# 
			# Please cite the following paper when you use the data generated by MemBrain:
			# 
			# Hongbin Shen, James J. Chou, MemBrain: Improving the Accuracy of Predicting
			#  Transmembrane Helices,  PLoS ONE, 2008, 6: e2399.
			# 
			# The query sequence:
			# MWELRSASFWRAIFAEFFATLFYVFFGLGASLRWAPGPLHVLQVALAFGL
			# ALATLVQAVGHISGAHVNPAVTFAFLVGSQMSLLRAICYVVAQLLGAVAG
			# AAVLYSVTPPAVRGNLALNTLHPGVSVGQATIVEIFLTLQFVLCIFATYD
			# ERRNGRLGSVALAVGFSLTLGHLFGMYYTGAGMNPARSFAPAILTRNFTN
			# HWVYWVGPVIGAGLGSLLYDFLLFPRLKSVSERLSILKGTRPSESNGQPE
			# VTGEPVELKTQAL
			# //
			# 
			# N-terminal signal peptide: NOT detected (based on user's knowledge)
			# //
			# 
			# Predicted Transmembrane Helix (TMH):
			# #1 TMH: 11-32
			# #2 TMH: 40-61
			# #3 TMH: 68-79 (possible half-TMH)
			# #4 TMH: 83-109
			# #5 TMH: 129-148
			# #6 TMH: 158-179
			# #7 TMH: 201-223
			# //
			# 
			# Propensity of residues forming TMH in the whole sequence:
			# -------------------------------------------------------------------------------------------------------------------
			# NOTE 1: Data is normalized to [0,1] and those larger than 0.4 are possible for forming TMH
			# NOTE 2: The TMH propensity plot of .gif format can be downloaded at http://www.csbio.sjtu.edu.cn/temp/11158/propensity.gif
			# NOTE 3: The TMH propensity plot of .eps format can be downloaded at http://www.csbio.sjtu.edu.cn/temp/11158/propensity.eps
			# -------------------------------------------------------------------------------------------------------------------
			# 1      0.000000      M
			# 2      0.000000      W
			# 3      0.000000      E
			# 4      0.000000      L

			if ($response2->is_success) {
				my $http_result2 = $response2->content;
				my @lines2 = split(/\n/, $http_result2);

				my $in_tmh_results = 0;
				my $seen_tmh_results = 0;
				my $tmh_string = '';

				foreach my $line2 (@lines2) {

					if ($seen_tmh_results == 0) {
						my $find_tmh_heading = index( $line2, 'Predicted Transmembrane Helix (TMH)' );
						if ($find_tmh_heading > -1) {
							$in_tmh_results = 1;
							$seen_tmh_results = 1;
						}
					} elsif ($in_tmh_results == 1) {
						my $find_tmh_end = index( $line2, '//' );
						if ($find_tmh_end > -1) {
							$in_tmh_results = 0;
						} else {
							my @bits = split( /\s/, $line2 );
							my $this_tmh = $bits[2];
							$this_tmh = trim($this_tmh);
							$tmh_string .= $this_tmh . ',';
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



sub trim {
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}

