#!/usr/bin/perl   -w

# perl perl_memsatsvm_post_and_save_link.pl -infile test1.fasta
# perl perl_memsatsvm_post_and_save_link.pl -infile test1.fasta -rmvggg

# If the -rmvggg flag is present, then remove leading and trailing GGG, HHH and AAA before submitting to prediction tool.

# will produce output file xxx.memsatsvm_link (eg. scrap.fasta.tmpred)
# will produce output file xxx.memsatsvm_html (eg. scrap.fasta.tmpred_html)
# will produce output file xxx.memsatsvm_redo (eg. scrap.fasta.tmpred_redo)

# The input to this program is a file containing multiple fasta sequences.
# Each fasta sequence is submitted to the MEMSAT-SVM server at http://bioinf.cs.ucl.ac.uk/psipred/?program=svmmemsat
# The REDO file looks just like the input file, and contains the sequences that always got rc 500 timeout,

# The TMpred program makes a prediction of membrane-spanning regions and their orientation. 
# The algorithm is based on the statistical analysis of TMbase, 
# a database of naturally occuring transmembrane proteins. 
# The prediction is made using a combination of several weight-matrices for scoring.



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
my $url_for_submit = 'http://bioinf.cs.ucl.ac.uk/psipred/submit';

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

my $output_file_link = "$fasta_file.memsatsvm_link";
my $output_file_2 = "$fasta_file.memsatsvm_html";
my $redo_file = "$fasta_file.memsatsvm_redo";

open( OUTFILE_LINK, ">$output_file_link") or
	die "Cannot open $output_file_link for writing\n";
if ($create_outfile2 == 1) {
	open( OUTFILE2, ">$output_file_2") or
		die "Cannot open $output_file_2 for writing\n";
}
open( REDOFILE, ">$redo_file") or
	die "Cannot open $redo_file for writing\n";

my $in  = Bio::SeqIO->new(-file => $fasta_file , '-format' => 'Fasta');

while ( my $seq = $in->next_seq() ) {
	process_1_fasta_seq($seq);
	sleep(10);
}

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
	if ($create_outfile2 == 1) {
		print OUTFILE2 "$time_str Processing $input_seqid\n\n";
	}

	my $ua = LWP::UserAgent->new;
	my $newagent = 'test search program (' . $ua->agent . ')';
	$ua->agent($newagent);

	my $input_fastaid_aaseq = "$input_fastaid\n$input_aaseq";

	my $request = POST ( $url_for_submit,
		Content => [ 	program => 'svmmemsat',
				sequence => $seq_and_index->{'seq'},
				complex => 'true',
				subject => $input_seqid
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

		# <h1 align="center">Your PSIPRED job has been submitted</h1>
		# <p align="center">Your job is in the queue under the name: 2zt9_G with the job ID: 378227 <p>
		# <p align="center"><img src="http://cms.cs.ucl.ac.uk/typo3/fileadmin/bioinf/templates/images/clock_2.gif" alt="clock"></p>
		# <p align="center">Jobs typically take around 5 minutes. If there is no response after 2 hours please email.</p>
		# <p align="center">Your results will be returned to this window once your job has completed. If you close the window or choose to navigate away you
		# can check the progress at <a href="http://bioinf.cs.ucl.ac.uk:80/psipred/result/378227">http://bioinf.cs.ucl.ac.uk:80/psipred/result/378227</a>.
		# Your results will be emailed to you once the job has completed</p>
		# 
		# <script type="text/javascript"><!--
		# setTimeout('Redirect()',20000);
		# function Redirect()
		# {
		#   location.href = 'http://bioinf.cs.ucl.ac.uk:80/psipred/result/378227';
		# }
		# // --></script>

		my $content = $response->content;

		my @lines = split(/\n/, $content);

		my $output_info = '';

		my $results_link = '';

		foreach my $line (@lines) {

			if (index($line, "location\.href = '") > -1) {
				my @bits = split( /location.href = '/, $line );
				my $bit1 = $bits[1];
				my @bits2 = split( /';/, $bit1 );
				$results_link = $bits2[0];
			}
		}

		if ($results_link eq '') {
			$results_link = '-';
		}
		my $output_line = $input_seqid . ':::::' . $results_link . ':::::';

		print OUTFILE_LINK "$output_line\n";
		close OUTFILE_LINK;
		open( OUTFILE_LINK, ">>$output_file_link") or
			die "Cannot open $output_file_link for writing\n";

		$results_link = '-';
		if ($results_link ne '-') {

			my $request2 = HTTP::Request->new( GET => $results_link );
			my $response2 = $ua->request($request2);

			#     <tr><td><h1>MEMSAT-SVM Prediction</h1></td></tr>
			#     <tr class="form_subtitle"><td><h4>Summary of MEMSAT-SVM Topology Analysis</h4></td></tr>
			# 
			# 
			# 
			# <tr><td>
			# <table>
			# 
			# 	  <tbody><tr>
			# 	  
			# 	    <td style="width: 20%;">Signal peptide</td>
			# 	    
			# 	    <td style="width: 20%;"><b>Not detected.</b></td>
			# 	    
			# 	  </tr>
			# 
			# 	  <tr>
			# 	  
			# 	    <td style="width: 20%;">Signal score</td>
			# 	    
			# 	    <td style="width: 20%;"><b>0</b></td>
			# 	    
			# 	  </tr>
			# 		
			# 	  <tr>
			# 	  
			# 	    <td style="width: 20%;">Topology</td>
			# 	    
			# 	    <td style="width: 20%;"><b>6-26</b></td>
			# 	    
			# 	  </tr>
			# 
			# 	  <tr>
			# 	  
			# 	    <td style="width: 20%;">Re-entrant helices</td>
			# 	    
			# 	    <td style="width: 20%;"><b>Not detected.</b></td>
			# 	    
			# 	  </tr>
			# 
			# 	  <tr>
			# 	  
			# 	    <td style="width: 20%;">Helix count</td>
			# 	    
			# 	    <td style="width: 20%;"><b>1</b></td>
			# 	    
			# 	  </tr>
			# 
			# 	  <tr>
			# 	  
			# 	    <td style="width: 20%;">N-terminal</td>
			# 	    
			# 	    <td style="width: 20%;"><b>out</b></td>
			# 	    
			# 	  </tr>
			# 
			# 	  <tr>
			# 	  
			# 	    <td style="width: 20%;">Score</td>
			# 	    
			# 	    <td style="width: 20%;"><b>2.47458</b></td>
			# 	    
			# 	  </tr>
			# 
			# </tbody></table>
			# </td></tr>
			# <tr><td><h3>Click to download these <a href="/bio_serf/getresultattached/7802963">results in plain text format</a></h3></td></tr>
			#      
			#   
			# 	<tr class="form_subtitle"><td><h4>MEMSAT-SVM Cartoon</h4></td></tr>
			# 	<tr><td><a href="/bio_serf/getresultattached/7802961"><img style="border: 1px none;" src="/bio_serf/getresultattached/7802961"></a></td></tr>
			#    
			#     <tr><td><h1>MEMSAT3 Prediction</h1></td></tr>
			#     <tr class="form_subtitle"><td><h4>Summary of MEMSAT3 Topology Analysis</h4></td></tr>
			# 
			# 
			#     <tr><td>
			#     <table style="width: 40%;">
			#       <tbody><tr><td><h4>Number</h4></td><td><h4>Type</h4></td><td><h4>Direction</h4></td><td><h4>Score</h4></td></tr>
			#       
			# 	  <tr style="padding: 3px; background: none repeat scroll 0% 0% rgb(204, 204, 204); text-align: center;"><td>1</td><td>helix</td><td>+</td><td>1.271</td></tr>
			# 
			# 	  <tr style="padding: 3px; background: none repeat scroll 0% 0% rgb(255, 255, 255); text-align: center;"><td>1</td><td>helix</td><td>-</td><td>31.901</td></tr>
			# 
			#     </tbody></table>
			#     </td></tr>
			#   
			#       <tr class="form_subtitle"><td><h4>Final MEMSAT3 Prediction</h4></td></tr>
			#       
			# 	  <tr><td>
			# 	  <table style="width: 40%;">
			# 	    <tbody><tr><td><h4>Segment</h4></td><td><h4>Range</h4></td><td><h4>Score</h4></td></tr>
			# 	  
			# 	    <tr style="padding: 3px; background: none repeat scroll 0% 0% rgb(204, 204, 204); text-align: center;"><td>1</td><td>(out) 7-26</td>
			# 	    
			# 	    <td>21.59</td></tr>
			# 	   
			# 	  </tbody></table>
			# 	  </td></tr>
			# 
			# 	  <tr class="form_subtitle"><td><h4>Diagram</h4></td></tr>
			# 	  <tr><td>
			# 	  <h3>Key:</h3>
			# 	  <table>
			# 	  <tbody><tr><td class="helix">X</td><td>Central transmembrane helix segment</td><td class="strand">S</td><td>Possible N-terminal signal peptide</td></tr>
			# 	  <tr><td class="inner_loop">+</td><td>Inside loop</td><td class="inner_cap">I</td><td>Inside helix cap</td></tr>
			# 	  <tr><td class="outer_loop">-</td><td>Outside loop</td><td class="outer_cap">O</td><td>Outside helix cap</td></tr>
			# 	  </tbody></table>
			# 	  </td></tr>
			# 
			# 	  
			# 	  <tr><td>
			# 	  <table>
			# 	  
			# 		<tbody><tr>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment">1</td>
			# 	
			# 		  <td class="alignment">0</td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment">2</td>
			# 	
			# 		  <td class="alignment">0</td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment">3</td>
			# 	
			# 		  <td class="alignment">0</td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		  <td class="alignment"> </td>
			# 	
			# 		</tr>
			# 	      
			# 
			# 				      
			# 		<tr>
			# 	
			# 		    <td class="outer_loop">-</td>
			# 		  
			# 		    <td class="outer_loop">-</td>
			# 		  
			# 		    <td class="outer_loop">-</td>
			# 		  
			# 		    <td class="outer_loop">-</td>
			# 		  
			# 		    <td class="outer_loop">-</td>
			# 		  
			# 		    <td class="outer_loop">-</td>
			# 		  
			# 		    <td class="outer_cap">O</td>
			# 		  
			# 		    <td class="outer_cap">O</td>
			# 		  
			# 		    <td class="outer_cap">O</td>
			# 		  
			# 		    <td class="outer_cap">O</td>
			# 		  
			# 		    <td class="outer_cap">O</td>
			# 		  
			# 		    <td class="helix">X</td>
			# 		  
			# 		    <td class="helix">X</td>
			# 		  
			# 		    <td class="helix">X</td>
			# 		  
			# 		    <td class="helix">X</td>
			# 		  
			# 		    <td class="helix">X</td>
			# 		  
			# 		    <td class="helix">X</td>
			# 		  
			# 		    <td class="helix">X</td>
			# 		  
			# 		    <td class="helix">X</td>
			# 		  
			# 		    <td class="helix">X</td>
			# 		  
			# 		    <td class="helix">X</td>
			# 		  
			# 		    <td class="inner_cap">I</td>
			# 		  
			# 		    <td class="inner_cap">I</td>
			# 		  
			# 		    <td class="inner_cap">I</td>
			# 		  
			# 		    <td class="inner_cap">I</td>
			# 		  
			# 		    <td class="inner_cap">I</td>
			# 		  
			# 		    <td class="inner_loop">+</td>
			# 		  
			# 		    <td class="inner_loop">+</td>
			# 		  
			# 		    <td class="inner_loop">+</td>
			# 		  
			# 		    <td class="inner_loop">+</td>
			# 		  
			# 		    <td class="inner_loop">+</td>
			# 		  
			# 		    <td class="inner_loop">+</td>
			# 		  
			# 		    <td class="inner_loop">+</td>
			# 		  
			# 		    <td class="inner_loop">+</td>
			# 		  
			# 		    <td class="inner_loop">+</td>
			# 		  
			# 		    <td class="inner_loop">+</td>
			# 		  
			# 		    <td class="inner_loop">+</td>
			# 		  
			# 		</tr>
			# 	      
			# 
			# 				      
			# 
			# 				      
			# 		<tr>
			# 	
			# 		  <td class="alignment">M</td>
			# 	       
			# 		  <td class="alignment">V</td>
			# 	       
			# 		  <td class="alignment">E</td>
			# 	       
			# 		  <td class="alignment">P</td>
			# 	       
			# 		  <td class="alignment">L</td>
			# 	       
			# 		  <td class="alignment">L</td>
			# 	       
			# 		  <td class="alignment">S</td>
			# 	       
			# 		  <td class="alignment">G</td>
			# 	       
			# 		  <td class="alignment">I</td>
			# 	       
			# 		  <td class="alignment">V</td>
			# 	       
			# 		  <td class="alignment">L</td>
			# 	       
			# 		  <td class="alignment">G</td>
			# 	       
			# 		  <td class="alignment">L</td>
			# 	       
			# 		  <td class="alignment">I</td>
			# 	       
			# 		  <td class="alignment">V</td>
			# 	       
			# 		  <td class="alignment">V</td>
			# 	       
			# 		  <td class="alignment">T</td>
			# 	       
			# 		  <td class="alignment">L</td>
			# 	       
			# 		  <td class="alignment">A</td>
			# 	       
			# 		  <td class="alignment">G</td>
			# 	       
			# 		  <td class="alignment">L</td>
			# 	       
			# 		  <td class="alignment">F</td>
			# 	       
			# 		  <td class="alignment">Y</td>
			# 	       
			# 		  <td class="alignment">A</td>
			# 	       
			# 		  <td class="alignment">A</td>
			# 	       
			# 		  <td class="alignment">Y</td>
			# 	       
			# 		  <td class="alignment">K</td>
			# 	       
			# 		  <td class="alignment">Q</td>
			# 	       
			# 		  <td class="alignment">Y</td>
			# 	       
			# 		  <td class="alignment">K</td>
			# 	       
			# 		  <td class="alignment">R</td>
			# 	       
			# 		  <td class="alignment">P</td>
			# 	       
			# 		  <td class="alignment">N</td>
			# 	       
			# 		  <td class="alignment">E</td>
			# 	       
			# 		  <td class="alignment">L</td>
			# 	       
			# 		  <td class="alignment">G</td>
			# 	       
			# 		  <td class="alignment">G</td>
			# 	       
			# 	       </tr>
			# 	      
			# 	     <tr><td colspan="35"><hr class="row_separator"></td></tr>
			# 	  
			# 	  </tbody></table>
			# 	  </td></tr>
			# 
			#   <tr><td>
			# <h3>Key:</h3>
			# 	  <table>
			# 	  <tbody><tr><td class="helix">X</td><td>Central transmembrane helix segment</td><td class="strand">S</td><td>Possible N-terminal signal peptide</td></tr>
			# 	  <tr><td class="inner_loop">+</td><td>Inside loop</td><td class="inner_cap">I</td><td>Inside helix cap</td></tr>
			# 	  <tr><td class="outer_loop">-</td><td>Outside loop</td><td class="outer_cap">O</td><td>Outside helix cap</td></tr>
			# 	  </tbody></table>
			#   </td></tr>
			# <tr><td><h3>Download these results in <a href="/bio_serf/getresultattached/7802964">plain text format</a></h3></td></tr>
			# 
			#   
			# 	<tr class="form_subtitle"><td><h4>MEMSAT Cartoon</h4></td></tr>
			# 	<tr><td><a href="/bio_serf/getresultattached/7802960"><img style="border: 1px none;" src="/bio_serf/getresultattached/7802960">

			if ($response2->is_success) {

				my $output_line = $input_seqid . ':::::';
				print OUTFILE "$output_line\n";

				my $http_result2 = $response2->content;
				my @lines2 = split(/\n/, $http_result2);

				my $need_to_output_current_result_lines = 0;

				foreach my $line2 (@lines2) {

					$line2 = trim($line2);
					if ($line2 ne '') {

						if ( (index($line2,"\<h1\>MEMSAT-SVM Prediction\<\/h1\>") > -1) || (index($line2,"\<h1\>MEMSAT3 Prediction\<\/h1\>") > -1) ) {
							$need_to_output_current_result_lines = 1;
						}
						if ( (index($line2,'Click to download these') > -1) || (index($line2,'Download these results in') > -1) ) {
							$need_to_output_current_result_lines = 1;
						}
						if ($need_to_output_current_result_lines == 1) {
							print OUTFILE "$line2\n";
						}
					}
				}

				close OUTFILE;
				open( OUTFILE, ">>$output_file_link") or
					die "Cannot open $output_file_link for writing\n";
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

