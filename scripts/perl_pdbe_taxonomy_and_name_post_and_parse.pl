#!/usr/bin/perl   -w

# perl perl_pdbe_taxonomy_and_name_post_and_parse.pl -infile opm_feb2012.unique_real_seq_opm_in_pdb.pdbid

# Read in a list of PDBIDs. 
# For each PDBID, read taxonomy from pdbe : http://www.ebi.ac.uk/pdbe-srv/view/entry/1zoy/taxonomy
# Output the taxonomy for each chain.
#
# Output the name for each chain, consisting of the crystal structure name plus the chain name

# input file :
# 1i78
# 1m56

# output file :
# 1i78_A:::::Species: Escherichia coli;;;;;Genus: Escherichia;;;;;Phylum: Proteobacteria;;;;;Superkingdom: Bacteria;;;;;cellular organisms;;;;;:::::
# 1i78_B:::::Species: Escherichia coli;;;;;Genus: Escherichia;;;;;Phylum: Proteobacteria;;;;;Superkingdom: Bacteria;;;;;cellular organisms;;;;;:::::
# 1zoy_A:::::Genus: Sus;;;;;Family: Suidae;;;;;Infraorder: Suina;;;;;Cetartiodactyla;;;;;Kingdom: Metazoa;;;;;Opisthokonta;;;;;Superkingdom: Eukaryota;;;;;cellular organisms;;;;;:::::
# 1zoy_B:::::Genus: Sus;;;;;Family: Suidae;;;;;Infraorder: Suina;;;;;Cetartiodactyla;;;;;Kingdom: Metazoa;;;;;Opisthokonta;;;;;Superkingdom: Eukaryota;;;;;cellular organisms;;;;;:::::

##################################################

use strict;
use Getopt::Long;
use LWP;
use HTTP::Request::Common;
use Data::Dumper;



##################################################

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
# http://www.ebi.ac.uk/pdbe-srv/view/entry/1zoy/taxonomy
my $url_1_for_submit = 'http://www.ebi.ac.uk/pdbe-srv/view/entry/';
my $url_2_for_submit = '/taxonomy';

# global variables

my $max_num_retries_of_query1 = 10;
my $count_input_pdbids = 0;
my $output_info;
my $referer = 'perl program';

# get the input parameters

my $input_file;
GetOptions( "infile=s" => \$input_file );
if (!defined $input_file) {
	die "Usage $0 -infile INPUTFILE\n";
}

# Open output files

my $output_file_taxonomy = "$input_file.pdbe_taxonomy";
my $output_file_name = "$input_file.pdbe_name";
my $output_file_skip = "$input_file.pdbe_skip";
my $output_file_2 = "$input_file.pdbe_taxonomy_html";
my $redo_file = "$input_file.pdbe_taxonomy_redo";

open( OUTFILE_TAX, ">$output_file_taxonomy") or
	die "Cannot open $output_file_taxonomy for writing\n";
open( OUTFILE_NAME, ">$output_file_name") or
	die "Cannot open $output_file_name for writing\n";
open( OUTFILE_SKIP, ">$output_file_skip") or
	die "Cannot open $output_file_skip for writing\n";
if ($create_outfile2 == 1) {
	open( OUTFILE2, ">$output_file_2") or
		die "Cannot open $output_file_2 for writing\n";
}
open( REDOFILE, ">$redo_file") or
	die "Cannot open $redo_file for writing\n";

# read input file

open INFILE, $input_file or die $!;
my @input_lines = <INFILE>;

# process each input line

foreach my $input_line (@input_lines) {
	chomp $input_line;
	$input_line = trim($input_line);
	if ($input_line ne '') {
		my @bits = split(/:::::/, $input_line);
		my $pdbid = $bits[0];
		$pdbid = substr( $pdbid, 0, 4 );
		process_1_pdbid($pdbid);
	}
}

# close the output files

#print OUTFILE "\n";
close OUTFILE_TAX;
close OUTFILE_NAME;
close OUTFILE_SKIP;
if ($create_outfile2 == 1) {
	print OUTFILE2 "\n";
	close OUTFILE2;
}
close REDOFILE;


##################################################
sub process_1_pdbid {

	my $pdbid = shift;

	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	$year += 1900;
	$mon += 1;
	$mon = sprintf("%02d", $mon);
	$mday = sprintf("%02d", $mday);
	$hour = sprintf("%02d", $hour);
	$min = sprintf("%02d", $min);
	$sec = sprintf("%02d", $sec);
	$time_str = "$year-$mon-$mday.$hour:$min:$sec";

	$count_input_pdbids++;
	print         "$time_str processing input number $count_input_pdbids, sequence $pdbid...\n";
	if ($create_outfile2 == 1) {
		print OUTFILE2 "$time_str Processing $pdbid\n\n";
	}

	my $ua = LWP::UserAgent->new;
	my $newagent = 'test search program (' . $ua->agent . ')';
	$ua->agent($newagent);

	my $url_for_submit = $url_1_for_submit . $pdbid . $url_2_for_submit;

	my $request = GET ( $url_for_submit,
		#Content => [ 	SEQUENCE => $form_aaseq
		#		],
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

	my $all_ok = 0;
	my $all_not_ok_msg = '';

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


		# http://www.ebi.ac.uk/pdbe-srv/view/entry/1zoy/taxonomy

#          <td class="contentsarea" id="contentsarea">
#                                                              
#                <h1>PDB entry 1zoy  
#                 
#                 </h1>                                                      
#                <span style="font-weight:bold">Crystal Structure of Mitochondrial Respiratory Complex II from porcine heart at 2.4 Angstroms</span>
#
#                <div style="height:15px;"></div>
#
#                
#                
#                <table class="content_table">
#                    <tr>
#                       <td class="content_header">Taxonomy</td>
#                    </tr>
#                    <tr>
#                      <td>
#
#
#                          <h2>
#                              FAD-binding protein  - Chain A</h2>
#                          <h3><a href="http://www.ebi.ac.uk/pdbe-srv/view/help#entity_src_nat" onmouseover="ajax_showTooltip('/pdbe-srv/view/help/entity_src_nat',this);return false;" onmouseout="ajax_hideTooltip()">Organism:</a></h3>
#                          <ol>
#                            
#                             <li>
#                                                             
#                            Sus scrofa 
#                                (<b>pig</b>)
#                              (<a target="external" href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&amp;id=9823">9823</a><a href="http://www.ebi.ac.uk/pdbe-srv/view/search/index/?Organism_tax_id=9823"><img src="www.ebi.ac.uk__pdbe-srv__view__entry__1zoy__taxonomy___files/srch.png" alt="search" border="0"></a>)
#                              <ol>
#                              
#                                <li style="font-size:95%">Genus: Sus 
#                                    
#                                    
#                                    (<a target="external" href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&amp;id=9822">9822</a><a href="http://www.ebi.ac.uk/pdbe-srv/view/search/index/?Hierarchy_tax_id=9822"><img src="www.ebi.ac.uk__pdbe-srv__view__entry__1zoy__taxonomy___files/srch.png" alt="search" border="0"></a>)
#                                     </li>
#                              
#                                <li style="font-size:95%">Family: Suidae 
#                                    - <b>pigs</b>
#
#                                <li style="font-size:95%">Opisthokonta 
#                                    
#                                    
#                                    (<a target="external" href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&amp;id=33154">33154</a><a href="http://www.ebi.ac.uk/pdbe-srv/view/search/index/?Hierarchy_tax_id=33154"><img src="www.ebi.ac.uk__pdbe-srv__view__entry__1zoy__taxonomy___files/srch.png" alt="search" border="0"></a>)
#                                     </li>
#                              
#                                <li style="font-size:95%">Superkingdom: Eukaryota 
#                                    - <b>eucaryotes</b>
#                                    
#                                    (<a target="external" href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&amp;id=2759">2759</a><a href="http://www.ebi.ac.uk/pdbe-srv/view/search/index/?Hierarchy_tax_id=2759"><img src="www.ebi.ac.uk__pdbe-srv__view__entry__1zoy__taxonomy___files/srch.png" alt="search" border="0"></a>)
#                                     </li>
#                              
#                                <li style="font-size:95%">cellular organisms 
#                                    
#                                    
#                                    (<a target="external" href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&amp;id=131567">131567</a><a href="http://www.ebi.ac.uk/pdbe-srv/view/search/index/?Hierarchy_tax_id=131567"><img src="www.ebi.ac.uk__pdbe-srv__view__entry__1zoy__taxonomy___files/srch.png" alt="search" border="0"></a>)
#                                     </li>
#                              
#                            </ol>
#                           
#
#                             </li>
#                            
#                          </ol>
#                        
#                        
#                        
#                          <h2>
#                              Iron-sulfur protein  - Chain B</h2>
#                          <h3><a href="http://www.ebi.ac.uk/pdbe-srv/view/help#entity_src_nat" onmouseover="ajax_showTooltip('/pdbe-srv/view/help/entity_src_nat',this);return false;" onmouseout="ajax_hideTooltip()">Organism:</a></h3>
#                          <ol>
#                            
#                             <li>
#                                                             
#                            Sus scrofa 
#                                (<b>pig</b>)

		# http://www.ebi.ac.uk/pdbe-srv/view/entry/1i78/taxonomy

#          <td class="contentsarea" id="contentsarea">
#                                                              
#                <h1>PDB entry 1i78  
#                 
#                 </h1>                                                      
#                <span style="font-weight:bold">CRYSTAL STRUCTURE OF OUTER MEMBRANE PROTEASE OMPT FROM ESCHERICHIA COLI</span>
#
#                <div style="height:15px;"></div>
#
#                
#                
#                <table class="content_table">
#                    <tr>
#                       <td class="content_header">Taxonomy</td>
#                    </tr>
#                    <tr>
#                      <td>
#
#
#                          <h2>
#                              PROTEASE VII  - Chains A, B</h2>
#                          <h3>Source organism</h3>
#                          <ol>
#                            
#                             <li>
#                                                             
#                            Escherichia coli K-12 
#                                
#                              (<a target="external" href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&amp;id=83333">83333</a><a href="http://www.ebi.ac.uk/pdbe-srv/view/search/index/?Organism_tax_id=83333"><img src="www.ebi.ac.uk__pdbe-srv__view__entry__1i78__taxonomy___files/srch.png" alt="search" border="0"></a>)
#                              <ol>
#                              
#                                <li style="font-size:95%">Species: Escherichia coli 
#                                    
#                                    
#                                    (<a target="external" href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&amp;id=562">562</a><a href="http://www.ebi.ac.uk/pdbe-srv/view/search/index/?Hierarchy_tax_id=562"><img src="www.ebi.ac.uk__pdbe-srv__view__entry__1i78__taxonomy___files/srch.png" alt="search" border="0"></a>)
#                                     </li>
#                              
#                                <li style="font-size:95%">Genus: Escherichia 
#                                    
#                                <li style="font-size:95%">Superkingdom: Bacteria 
#                                    - <b>eubacteria</b>
#                                    
#                                    (<a target="external" href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&amp;id=2">2</a><a href="http://www.ebi.ac.uk/pdbe-srv/view/search/index/?Hierarchy_tax_id=2"><img src="www.ebi.ac.uk__pdbe-srv__view__entry__1i78__taxonomy___files/srch.png" alt="search" border="0"></a>)
#                                     </li>
#                              
#                                <li style="font-size:95%">cellular organisms 



		my @lines = split(/\n/, $content);

		my $pdbid_name = '';
		my $taxonomy_output = '';
		my @chains_id;
		my @chains_name;
		my $in_a_chain = 0;
		my $pdbe_didnt_return_much_output = 0;

		foreach my $line (@lines) {
			chomp $line;
			$line = trim($line);

			#print "$line\n\n\n\n";

			if ($line =~ /All the polymer chains in this entry are synthetic./) {
				$pdbe_didnt_return_much_output = 1;
				$in_a_chain = 0;
				my $output_line = "$pdbid : All the polymer chains in this entry are synthetic.\n";
				print OUTFILE_SKIP $output_line;
				close OUTFILE_SKIP;
				open( OUTFILE_SKIP, ">>$output_file_skip") or
					die "Cannot open $output_file_skip for appending\n";
			}

			if ($line =~ /The experimental method for this entry is THEORETICAL MODEL./) {
				$pdbe_didnt_return_much_output = 1;
				$in_a_chain = 0;
				my $output_line = "$pdbid : The experimental method for this entry is THEORETICAL MODEL.\n";
				print OUTFILE_SKIP $output_line;
				close OUTFILE_SKIP;
				open( OUTFILE_SKIP, ">>$output_file_skip") or
					die "Cannot open $output_file_skip for appending\n";
			}

			if (($line ne '') and ($pdbe_didnt_return_much_output == 0)) {

				if ($line =~ /<span style="font-weight:bold">/) {
					# <span style="font-weight:bold">Crystal Structure of Mitochondrial Respiratory Complex II from porcine heart at 2.4 Angstroms</span>
					my @bits1 = split(/<span style="font-weight:bold">/, $line);
					my $bit2 = trim($bits1[1]);
					my @bits3 = split(/<\/span>/, $bit2);
					$pdbid_name = trim($bits3[0]);

				} elsif ($line =~ /Chain/) {

					# if already in a chain, terminate it
					if ($in_a_chain == 1) {
						my $i = 0;
						foreach my $this_chain (@chains_id) {
							my $this_chain_name = $chains_name[$i];
							if ($pdbid_name eq '') {
								$pdbid_name = '-';
							}
							if ($taxonomy_output eq '') {
								$taxonomy_output = '-';
							}
							my $output_line = $pdbid . '_' . $this_chain . ':::::' . $pdbid_name . ':::::' . $this_chain_name . ':::::' . "\n";
							print OUTFILE_NAME $output_line;
							$output_line = $pdbid . '_' . $this_chain . ':::::' . $taxonomy_output . ':::::' . "\n";
							print OUTFILE_TAX $output_line;
							$i++;
						}
					}

					# start a new chain
					$in_a_chain = 1;
					$taxonomy_output = '';
					@chains_id = ();
					@chains_name = ();

					if ($line =~ /Chains/) {
						# PROTEASE VII  - Chains A, B</h2>
						my @bits1 = split(/<\/h2>/, $line);
						my $bit2 = trim($bits1[0]);
						my @bits3 = split(/- Chains /, $bit2);
						my $this_chain_name = trim($bits3[0]);
						my $bit4 = $bits3[1];
						my @bits5 = split(/, /, $bit4);
						foreach my $this_chain (@bits5) {
							push( @chains_id, $this_chain );
							push( @chains_name, $this_chain_name );
						}
					} else {
						# FAD-binding protein  - Chain A</h2>
						my @bits1 = split(/<\/h2>/, $line);
						my $bit2 = trim($bits1[0]);
						my @bits3 = split(/- Chain /, $bit2);
						my $this_chain_name = trim($bits3[0]);
						my $bit4 = $bits3[1];
						push( @chains_id, $bit4 );
						push( @chains_name, $this_chain_name );
					}

				} elsif ($line =~ /<li style="font-size:95%">/) {
					# <li style="font-size:95%">Genus: Sus
					my @bits1 = split(/<li style="font-size:95%">/, $line);
					my $bit2 = trim($bits1[1]);
					$taxonomy_output .= $bit2 . ';;;;;';
				}
			}
		}

		# if in a chain, terminate it when read end of html input file
		if ($in_a_chain == 1) {
			my $i = 0;
			foreach my $this_chain (@chains_id) {
				my $this_chain_name = $chains_name[$i];
				if ($pdbid_name eq '') {
					$pdbid_name = '-';
				}
				if ($taxonomy_output eq '') {
					$taxonomy_output = '-';
				}
				my $output_line = $pdbid . '_' . $this_chain . ':::::' . $pdbid_name . ':::::' . $this_chain_name . ':::::' . "\n";
				print OUTFILE_NAME $output_line;
				$output_line = $pdbid . '_' . $this_chain . ':::::' . $taxonomy_output . ':::::' . "\n";
				print OUTFILE_TAX $output_line;
				$i++;
			}
		}


	} else {
		print "     Had $max_num_retries_of_query1 timeouts - got no results.\n";
		$all_not_ok_msg = "Had $max_num_retries_of_query1 timeouts - got no results.";
	}

	if ($all_ok == 0) {
		print REDOFILE "\>$pdbid\n\n";
	}

	close OUTFILE_TAX;
	close OUTFILE_NAME;
	open( OUTFILE_TAX, ">>$output_file_taxonomy") or
		die "Cannot open $output_file_taxonomy for appending\n";
	open( OUTFILE_NAME, ">>$output_file_name") or
		die "Cannot open $output_file_name for appending\n";
}

##################################################
# Perl trim function to remove whitespace from the start and end of the string
sub trim {
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}

