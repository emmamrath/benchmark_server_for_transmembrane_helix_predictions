#!/usr/bin/perl   -w

# perl perl_memsatsvm_get_link_and_save_results.pl -infile alpha_polytopic_bitopic_noblanks.fasta.memsatsvm1.memsatsvm_link

# will produce output file xxx.memsatsvm_results
# will produce output file xxx.memsatsvm_html
# will produce output file xxx.memsatsvm_redo

# The input to this program is a file containing multiple fasta sequences.
# Each fasta sequence was submitted to the MEMSAT-SVM server at http://bioinf.cs.ucl.ac.uk/psipred/?program=svmmemsat
# and the HTML results gave a link that was saved into the input file of this program.
# This program will do an HTTP GET of that link, to get the actual results.
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

my $input_file;
GetOptions( "infile=s" => \$input_file );
if (!defined $input_file) {
	die "Usage $0 -infile 1LINEfile\n";
}

# Open output files

my $output_file = "$input_file.memsatsvm_results";
my $output_file_2 = "$input_file.memsatsvm_html2";
my $redo_file = "$input_file.memsatsvm_redo2";

open( OUTFILE, ">$output_file") or
	die "Cannot open $output_file for writing\n";
if ($create_outfile2 == 1) {
	open( OUTFILE2, ">$output_file_2") or
		die "Cannot open $output_file for writing\n";
}
open( REDOFILE, ">$redo_file") or
	die "Cannot open $output_file for writing\n";

# processing

open INFILE, $input_file or die $!;
my @input_lines = <INFILE>;

foreach my $input_line (@input_lines) {

	chomp( $input_line );
	$input_line = trim($input_line);
	if ($input_line ne '') {

		# write the fasta-id and fasta-sequence to the temporary file,
		# so it can be input to the tmhmm command

		my @bits = split(/:::::/, $input_line);

		my $this_id = $bits[0];
		my $this_link = $bits[1];

		if ($this_link ne '-') {

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
			print         "$time_str processing input number $count_input_seqs, sequence $this_id...\n";
			# print OUTFILE "$time_str Processing $this_id\n";
			# print OUTFILE "$input_aaseq\n";
			if ($create_outfile2 == 1) {
				print OUTFILE2 "$time_str Processing $this_id\n\n";
			}

			my $ua = LWP::UserAgent->new;
			my $newagent = 'test search program (' . $ua->agent . ')';
			$ua->agent($newagent);

			my $request = HTTP::Request->new( GET => $this_link );
			my $response = $ua->request($request);

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

			if ($response->is_success) {

				my $output_line = $this_id . ':::::';
				print OUTFILE "$output_line\n";

				my $http_result2 = $response->content;
				my @lines2 = split(/\n/, $http_result2);

				my $need_to_output_current_result_lines = 0;

				foreach my $line2 (@lines2) {

					$line2 = trim($line2);
					if ($line2 ne '') {

						if ( (index($line2,"\<h1\>MEMSAT-SVM Prediction\<\/h1\>") > -1) || (index($line2,"\<h1\>MEMSAT3 Prediction\<\/h1\>") > -1) ) {
							$need_to_output_current_result_lines = 1;
						}
						if ( (index($line2,'Click to download these') > -1) || (index($line2,'Download these results in') > -1) ) {
							$need_to_output_current_result_lines = 0;
						}
						if ($need_to_output_current_result_lines == 1) {
							print OUTFILE "$line2\n";
						}
					}
				}

				close OUTFILE;
				open( OUTFILE, ">>$output_file") or
					die "Cannot open $output_file for writing\n";
			}
		}
	}
}

print OUTFILE "\n";
close OUTFILE;
if ($create_outfile2 == 1) {
	print OUTFILE2 "\n";
	close OUTFILE2;
}
close REDOFILE;



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

