#!/usr/bin/perl   -w

# perl perl_minnou_post_and_save_link.pl -infile identity_30_graded.1line -rmvggg
# perl perl_minnou_post_and_save_link.pl -infile alpha_polytopic_bitopic_noblanks.1line
# perl perl_minnou_post_and_save_link.pl -infile beta_barrel_noblanks.1line
# perl perl_minnou_post_and_save_link.pl -infile soluble.1line
# perl perl_minnou_post_and_save_link.pl -infile test1.1line

# If the -rmvggg flag is present, then remove leading and trailing GGG, HHH and AAA before submitting to prediction tool.

# will produce output file xxx.minnou_link (eg. identity_30_graded.1line.minnou_link)

# Sequences were submitted to the minnou server http://minnou.cchmc.org/
# to predict transmembrane helices.

# Example of infile :
# 1uaz_A:::::>1uaz_A:::::TAAVGADLLGDGRPETLWLGIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIGLTEVQVGSEMLDIYYARYADWLFTTPLLLLDLALLAKVDRVSIGTLVGVDALMIVTGLVGALSHTPLARYTWWLFSTICMIVVLYFLATSLRAAAKERGPEVASTFNTLTALVLVLWTAYPILWIIGTEGAGVVGLGIETLLFMVLDVTAKVGFGFILLRSRAILGDTEAPEPSAGAEASAAD:::::
# 3ddl_A:::::>3ddl_A:::::MLQELPTLTPGQYSLVFNMFSFTVATMTASFVFFVLARNNVAPKYRISMMVSALVVFIAGYHYFRITSSWEAAYALQNGMYQPTGELFNDAYRYVDWLLTVPLLTVELVLVMGLPKNERGPLAAKLGFLAALMIVLGYPGEVSENAALFGTRGLWGFLSTIPFVWILYILFTQLGDTIQRQSSRVSTLLGNARLLLLATWGFYPIAYMIPMAFPEAFPSNTPGTIVALQVGYTIADVLAKAGYGVLIYNIAKAKSEEEGFNVSEMVEPATASA:::::
# 1f6r_A:::::>1f6r_A:::::EQLTKCEVFRELKDLKGYGGVSLPEWVCTTFHTSGYDTQAIVQNNDSTEYGLFQINNKIWCKDDQNPHSSNICNISCDKFLDDDLTDDIMCVKKILDKVGINYWLAHKALCSEKLDQWLCEKL:::::

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
my $url_for_submit = 'http://minnou.cchmc.org/cgi-bin/predict.cgi';

# global variables

my $max_num_retries_of_query1 = 1;
my $count_input_seqs = 0;
my $referer = 'Client user agent';

my $input_file;
my $input_rmvggg;
GetOptions( "infile=s" => \$input_file, "rmvggg" => \$input_rmvggg );
if (!defined $input_file) {
	die "Usage $0 -infile INFILE\n";
}

# Open output files

my $output_file = "$input_file.minnou_link";
my $output_file_2 = "$input_file.minnou_html";
my $redo_file = "$input_file.minnou_redo";

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

my $jobid = 1633333920;

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

	# <form action="frame.php" method="POST" enctype="multipart/form-data">
	# <input TYPE="hidden" name="MAX_FILE_SIZE" value="100000000">
	# <br><input NAME="tmpsubmittedfile" TYPE="file">
	# <br><textarea NAME="sequences" ROWS=5 COLS=50>
	# </textarea><br><br><br><br><br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
	# <input TYPE="submit" VALUE="Submit job">&nbsp;&nbsp;
	# <input type="reset" value="Reset">
	# <input TYPE=hidden name="jobid" value=1416980517>
	# <input TYPE=hidden name=p value="main">

	$jobid = $jobid + 1;
	my $input_field = "$fastaid\n$aaseq";
	my $request = POST ( $url_for_submit,
		Content => [ 	seqName => $fastaid,
				aaSeq => $seq_and_index->{'seq'},
				RLApred => '',
				RLAmodel => '1',
				SABLEres => '',
				Email => ''
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
		my $output_line = $seqid . ':::::ERROR:::::';
		my $full_link = '';
		my $partial_link = '';

		foreach my $line (@lines) {

			#print "$line\n";
			print OUTFILE2 "$line\n";

			# To view results, <a href="http://minnou.cchmc.org/cgi-bin/view_results.cgi?ID=df7498ce559bdc0">click here</a>.

			my $search_string = 'To view results';
			my $search_string_result = index( $line, $search_string );
			if ( $search_string_result != -1 ) {
				my @bits = split( /To view results, \<a href=\"/, $line );
				my @bits2 = split( /\"\>click here\<\/a\>/, $bits[1] );
				$partial_link = $bits2[0];
				$full_link = $partial_link;

			}
		}

		if ($partial_link eq '') {
			$full_link = '-';
		}

		$output_line = $seqid . ':::::' . $full_link . ':::::' . $aaseq . ':::::';
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




