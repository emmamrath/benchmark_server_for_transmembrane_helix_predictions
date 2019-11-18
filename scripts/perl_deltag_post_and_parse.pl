#!/usr/bin/perl   -w

# perl perl_deltag_post_and_parse.pl -infile identity_30_graded.1line -rmvggg
# perl perl_deltag_post_and_parse.pl -infile scrap.fasta.1line

# If the -rmvggg flag is present, then remove leading and trailing GGG, HHH and AAA before submitting to prediction tool.

# will produce output file xxx.deltag (eg. identity_30_graded.1line.deltag)

# Sequences were submitted to the deltaG server http://dgpred.cbr.su.se/index.php?p=fullscan
# to predict transmembrane helices.

# Example of infile :
# 1uaz_A:::::>1uaz_A:::::TAAVGADLLGDGRPETLWLGIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIGLTEVQVGSEMLDIYYARYADWLFTTPLLLLDLALLAKVDRVSIGTLVGVDALMIVTGLVGALSHTPLARYTWWLFSTICMIVVLYFLATSLRAAAKERGPEVASTFNTLTALVLVLWTAYPILWIIGTEGAGVVGLGIETLLFMVLDVTAKVGFGFILLRSRAILGDTEAPEPSAGAEASAAD:::::
# 3ddl_A:::::>3ddl_A:::::MLQELPTLTPGQYSLVFNMFSFTVATMTASFVFFVLARNNVAPKYRISMMVSALVVFIAGYHYFRITSSWEAAYALQNGMYQPTGELFNDAYRYVDWLLTVPLLTVELVLVMGLPKNERGPLAAKLGFLAALMIVLGYPGEVSENAALFGTRGLWGFLSTIPFVWILYILFTQLGDTIQRQSSRVSTLLGNARLLLLATWGFYPIAYMIPMAFPEAFPSNTPGTIVALQVGYTIADVLAKAGYGVLIYNIAKAKSEEEGFNVSEMVEPATASA:::::

# http://www.cbr.su.se/DGpred/


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
my $url_for_submit = 'http://dgpred.cbr.su.se/results.php?program=fullscan';

# global variables

my $max_num_retries_of_query1 = 10;
my $count_input_seqs = 0;
my $referer = 'Client user agent';

my $input_file;
my $input_rmvggg;
GetOptions( "infile=s" => \$input_file, "rmvggg" => \$input_rmvggg );
if (!defined $input_file) {
	die "Usage $0 -infile INFILE\n";
}

# Open output files

my $output_file = "$input_file.deltag";
my $output_file_2 = "$input_file.deltag_html";
my $redo_file = "$input_file.deltag_redo";

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
	my $newagent = 'Client user agent (' . $ua->agent . ')';
	$ua->agent($newagent);

	# <form action="/results.php?program=fullscan" method="post" enctype="multipart/form-data">
	# <textarea name="seq" rows="10" cols="80"></textarea><br><br>
	# <a href="#" onClick="javascript:hideShow('options');return false;" style="text-decoration:none; color:#333333"><b>Options</b> <img id="arrow" src='arrow90.png' border=0></a>
	# <table cellspacing=0 cellpadding=0 border=0 id="options">
	# <tr><td colspan=2><img src='0.gif' width='1' height='10'></td></tr>
	# <tr><td>Helix min length: &nbsp;&nbsp;&nbsp;</td><td align="center"> <input type="text" name="Lmin" size="2" value="19"></input></td></tr>
	# <tr><td>Helix max length: &nbsp;&nbsp;&nbsp;</td><td align="center"> <input type="text" name="Lmax" size="2" value="23"></input></td></tr>
	# <tr><td>Length correction: &nbsp;&nbsp;&nbsp;</td><td align="center"> <input type="checkbox" name="with_length" checked></input></td></tr>
	# <tr><td>Show graphics: &nbsp;&nbsp;&nbsp;</td><td align="center"> <input type="checkbox" name="graphics" checked></input></td></tr>
	# </table>
	# <br><br>
	# <input type="submit" value="Submit">
	# <input type="reset" value="Clear">
	# </form>

	my $input_field = "$fastaid\n$aaseq";
	my $request = POST ( $url_for_submit,
		Content => [ 	seq => $seq_and_index->{'seq'},
				'Lmin' => '19',
				'Lmax' => '23',
				'with_length' => 'on'
				#'graphics' => '0'
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

		my $output_line = $seqid . ':::::-:::::';
		my $output_seq = '';
		my $seen_start = 0;
		my $output_1line = '';

		foreach my $line (@lines) {

			# <table cellspacing="0" cellpadding="0" border="0" width="720" class="TableGrey">
			# <tr><td align="left"><table cellspacing=0 cellpadding=5 border=0 class="TableWhite" width='100%'>
			# <tr><td colspan=3><big>&gt;3o7q_A</big></td></tr>
			# <tr><td colspan=2 valign="top"><img src="0.gif" height=1 width=200><br><table cellspacing=0 cellpadding=0 border=0><tr><td colspan=9><i><b>Predicted TM helices:</b><br></i></td></tr><tr><td align=center><i>Position</i></td><td>&nbsp;&nbsp;&nbsp;</td><td align=center><i>Length</i></td><td align=center>&nbsp;&nbsp;&nbsp;</td><td><i>Predicted &Delta;G</i></td><td>&nbsp;&nbsp;&nbsp;</td><td align=left><i>Sequence</i></td><td><img src='0.gif' height=1 width=100></td><td align=left>&nbsp;</td></tr><tr><td align=center>22-41</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>20</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>-1.518</td><td>&nbsp;</td><td>RSYIIPFALLCSLFFLWAVA</td><td>&nbsp;</td><td><a href="analyze.php?with_length=on&seq=RSYIIPFALLCSLFFLWAVA" target="_blank" style='text-decoration:none'  onClick="javascript:analyze('RSYIIPFALLCSLFFLWAVA');return false;"><span style='text-decoration: none;color: #800000;'>[analyze]</span></a></td></tr>
			# <tr><td align=center>69-88</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>20</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>0.909</td><td>&nbsp;</td><td>AFYFGYFIIPIPAGILMKKL</td><td>&nbsp;</td><td><a href="analyze.php?with_length=on&seq=AFYFGYFIIPIPAGILMKKL" target="_blank" style='text-decoration:none'  onClick="javascript:analyze('AFYFGYFIIPIPAGILMKKL');return false;"><span style='text-decoration: none;color: #800000;'>[analyze]</span></a></td></tr>
			# <tr><td align=center>90-112</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>23</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>-0.506</td><td>&nbsp;</td><td>YKAGIITGLFLYALGAALFWPAA</td><td>&nbsp;</td><td><a href="analyze.php?with_length=on&seq=YKAGIITGLFLYALGAALFWPAA" target="_blank" style='text-decoration:none'  onClick="javascript:analyze('YKAGIITGLFLYALGAALFWPAA');return false;"><span style='text-decoration: none;color: #800000;'>[analyze]</span></a></td></tr>
			# <tr><td align=center>114-134</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>21</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>-1.690</td><td>&nbsp;</td><td>IMNYTLFLVGLFIIAAGLGCL</td><td>&nbsp;</td><td><a href="analyze.php?with_length=on&seq=IMNYTLFLVGLFIIAAGLGCL" target="_blank" style='text-decoration:none'  onClick="javascript:analyze('IMNYTLFLVGLFIIAAGLGCL');return false;"><span style='text-decoration: none;color: #800000;'>[analyze]</span></a></td></tr>
			# <tr><td align=center>160-179</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>20</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>1.655</td><td>&nbsp;</td><td>TFNSFGAIIAVVFGQSLILS</td><td>&nbsp;</td><td><a href="analyze.php?with_length=on&seq=TFNSFGAIIAVVFGQSLILS" target="_blank" style='text-decoration:none'  onClick="javascript:analyze('TFNSFGAIIAVVFGQSLILS');return false;"><span style='text-decoration: none;color: #800000;'>[analyze]</span></a></td></tr>
			# <tr><td align=center>210-231</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>22</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>-3.553</td><td>&nbsp;</td><td>TPYMIIVAIVLLVALLIMLTKF</td><td>&nbsp;</td><td><a href="analyze.php?with_length=on&seq=TPYMIIVAIVLLVALLIMLTKF" target="_blank" style='text-decoration:none'  onClick="javascript:analyze('TPYMIIVAIVLLVALLIMLTKF');return false;"><span style='text-decoration: none;color: #800000;'>[analyze]</span></a></td></tr>
			# <tr><td align=center>260-282</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>23</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>1.406</td><td>&nbsp;</td><td>WRWAVLAQFCYVGAQTACWSYLI</td><td>&nbsp;</td><td><a href="analyze.php?with_length=on&seq=WRWAVLAQFCYVGAQTACWSYLI" target="_blank" style='text-decoration:none'  onClick="javascript:analyze('WRWAVLAQFCYVGAQTACWSYLI');return false;"><span style='text-decoration: none;color: #800000;'>[analyze]</span></a></td></tr>
			# <tr><td align=center>296-318</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>23</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>-0.012</td><td>&nbsp;</td><td>FAANYLTGTMVCFFIGRFTGTWL</td><td>&nbsp;</td><td><a href="analyze.php?with_length=on&seq=FAANYLTGTMVCFFIGRFTGTWL" target="_blank" style='text-decoration:none'  onClick="javascript:analyze('FAANYLTGTMVCFFIGRFTGTWL');return false;"><span style='text-decoration: none;color: #800000;'>[analyze]</span></a></td></tr>
			# <tr><td align=center>327-345</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>19</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>-1.267</td><td>&nbsp;</td><td>VLAAYALIAMALCLISAFA</td><td>&nbsp;</td><td><a href="analyze.php?with_length=on&seq=VLAAYALIAMALCLISAFA" target="_blank" style='text-decoration:none'  onClick="javascript:analyze('VLAAYALIAMALCLISAFA');return false;"><span style='text-decoration: none;color: #800000;'>[analyze]</span></a></td></tr>
			# <tr><td align=center>349-369</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>21</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>1.212</td><td>&nbsp;</td><td>VGLIALTLCSAFMSIQYPTIF</td><td>&nbsp;</td><td><a href="analyze.php?with_length=on&seq=VGLIALTLCSAFMSIQYPTIF" target="_blank" style='text-decoration:none'  onClick="javascript:analyze('VGLIALTLCSAFMSIQYPTIF');return false;"><span style='text-decoration: none;color: #800000;'>[analyze]</span></a></td></tr>
			# <tr><td align=center>382-404</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>23</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>1.485</td><td>&nbsp;</td><td>YGSSFIVMTIIGGGIVTPVMGFV</td><td>&nbsp;</td><td><a href="analyze.php?with_length=on&seq=YGSSFIVMTIIGGGIVTPVMGFV" target="_blank" style='text-decoration:none'  onClick="javascript:analyze('YGSSFIVMTIIGGGIVTPVMGFV');return false;"><span style='text-decoration: none;color: #800000;'>[analyze]</span></a></td></tr>
			# <tr><td align=center>413-432</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>20</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>-0.493</td><td>&nbsp;</td><td>TAELIPALCFAVIFIFARFR</td><td>&nbsp;</td><td><a href="analyze.php?with_length=on&seq=TAELIPALCFAVIFIFARFR" target="_blank" style='text-decoration:none'  onClick="javascript:analyze('TAELIPALCFAVIFIFARFR');return false;"><span style='text-decoration: none;color: #800000;'>[analyze]</span></a></td></tr>
			# </table><br></td>
			# </tr></table><br><br><br>
			# </td></tr></table>

			$line = trim($line);
			if ($line ne '') {

				if (index($line, 'Predicted TM helices:') > -1) {
					$seen_start = 1;
				}
				if ($seen_start == 1) {
					if (index($line, '</table><br></td>') > -1) {
						$seen_start = 0;
					} else {
						$output_1line .= $line;
					}
				}
			}
		}

		#print "\n$output_1line\n\n"; 

		my @lines2 = split(/\<tr\>\<td align=center\>/, $output_1line);

		foreach my $line (@lines2) {
			$line = trim($line);
			if ($line ne '') {

				if ((index($line, 'Predicted TM helices:') == -1) && (index($line, 'Position') == -1)) {

					my @bits = split( /\<\/td>/, $line );
					my $bit1 = $bits[0];
					$output_seq .= $bit1 . ',';
				}
			}
		}

		if ($output_seq eq '') {
			$output_seq = '-';
		}
		$output_line = $seqid . ':::::' . $output_seq . ':::::';
		print OUTFILE "$output_line\n";

		close OUTFILE;
		open( OUTFILE, ">>$output_file") or
			die "Cannot open $output_file for rewriting\n";

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




