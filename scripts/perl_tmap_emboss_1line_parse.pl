#!/usr/bin/perl   -w

# perl perl_tmap_emboss_1line_parse.pl -infile scrap.1line -rmvggg
# perl perl_tmap_emboss_1line_parse.pl -infile all.1line -rmvggg

# If the -rmvggg flag is present, then remove leading and trailing GGG, HHH and AAA before submitting to prediction tool.

# will produce output file xxx.tmap
# will also produce a temporary file xxx.tmhmm_temp that will be left there permanently (eg. scrap1_seqids.unique.1line.tmhmm_temp)
# will also produce a temporary file xxx.tmhmm_temp2 that will be left there permanently (eg. scrap1_seqids.unique.1line.tmhmm_temp2)

# The input to this program is a file containing seqids, fasta-ids, and aa-seqs, 1 line per seqid.
# The output has, for each seqid, the tmhmm results.
# The call to tmhmm needs a input file containing the 1 fasta sequence,
# which is why this program recreates the temporary file for each sequence processed.
# The second temporary file will contain the output returned by tmhmm.

# Example of input file:
# AAA26845:::::>gi|153550|gb|AAA26845.1| penicillin-binding protein (PBPB2):::::AVIASISKEMPGISISTSWDRKILETSLSSIVGSVSSEKAGLPAEEAETYLKKGYSLNDRVGTSYLEKQYEETLQGKRSVKEIHLDKYGNMESVENIEDGTKGNNIKLTIDLSFQDSVDALLKSYFNSELGNGGAKYSEGVYAVALNPKTGAVLSMSGIKHDLKTGELTPDSLGTVTNVFVPGSVVKAATISSGWENGVLSGNQTLTDQSIVFQGSAPINSWYTAFSRPMPITAVQALEYSSNAYMVQTALGLMGQTYQPNMFVGTSNLESAMGKLRSTFGEYGLGSATGIDLPDESTGFIPKEYSFANFITNAFGQFDNYTPMQLAQYVATIANDGVRVAPRIVEGIYGNNDKGGLGGLIQQLQPTEMNKVNISDSDMSVLHQGFYQVAHGTSGLTTGRAFSNGAAVSISGKTGTAESYVAGGQEANNTNAVAYAPSDNPQIAVAVVFPHNTNLTNGVGPSIARDIINLYNQHHPMN
# AAA26872:::::>gi|153613|gb|AAA26872.1| DpnII DNA methylase:::::MKIKEIKKVTLQPFTKWTGGKRQLLPVIRELIPKTYNRYFEPFVGGGALFFDLAPKDAVINDFNAELINCYQQIKDNPQELIEILKVHQEYNSKEYYLDLRSADRDERIDMMSEVQRAARILYMLRVNFNGLYRVNSKNQFNVPYGRYKNPKIVDEELISAISVYINNNQLEIKVGDFEKAIVDVRTGDFVYFDPPYIPLSETSAFTSYTHEGFSFADQVRLRDAFKRLSDTGAYVMLSNSSSALVEELYKDFNIHYVEATRTNGAKSSSRGKISEIIVTNYEK



use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use IO::File;
use XML::Parser;
use Data::Dumper;

# Check to see if -infile was specified

my $input_file;
my $input_rmvggg;
GetOptions( "infile=s" => \$input_file, "rmvggg" => \$input_rmvggg );
if (!defined $input_file) {
	die "Usage $0 -infile 1LINESfile\n";
}

# Open output file
my $output_file = "$input_file.tmap";
my $temp_file = "$input_file.tmap_temp";
my $temp_file_2 = "$input_file.tmap_temp2";
open( OUTFILE, ">$output_file") or
	die "Cannot open $output_file for writing\n";

# processing

open INFILE, $input_file or die $!;
my @input_lines = <INFILE>;

my $upto = 0;

foreach my $input_line (@input_lines) {

	chomp( $input_line );
	$input_line = trim($input_line);
	if ($input_line ne '') {

		# write the fasta-id and fasta-sequence to the temporary file,
		# so it can be input to the tmhmm command

		my @bits = split(/:::::/, $input_line);

		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		$year += 1900;
		$mon += 1;
		$mon = sprintf("%02d", $mon);
		$mday = sprintf("%02d", $mday);
		$hour = sprintf("%02d", $hour);
		$min = sprintf("%02d", $min);
		$sec = sprintf("%02d", $sec);
		my $time_str = "$year-$mon-$mday.$hour:$min:$sec";

		$upto++;
		my $seqid = $bits[0];
		print "$time_str : Processing sequence number $upto : $seqid\n";

		my $fasta_id = $bits[1];
		my $output_redo = '-';
		my $input_aaseq = $bits[2];

		my $input_aaseq_fix = $input_aaseq;
		$input_aaseq_fix =~ s/X/G/g; # convert invalid residue to valid, make sure program doesn't ignore residue and count residue positions incorrectly 

		my %seq_and_index_array;
		my $seq_and_index = \%seq_and_index_array;
		$seq_and_index->{'seq'} = $input_aaseq_fix;
		$seq_and_index->{'index'} = 0;
		if (defined($input_rmvggg)) {
			$seq_and_index = find_leading_trailing_GGG_HHH_AAA( $input_aaseq_fix );
		}

		my $input_msf = '';
		my $input_msf_name = 'XXXX_YYYYY';
		my $input_msf_length = length( $seq_and_index->{'seq'} );
		

		my $aa_seq = '';
		my $max_str_len = 50;
		my $remaining_seq = $seq_and_index->{'seq'};
		my $keep_looping = 1;
		while ($keep_looping == 1) {

			my $str_len = length $remaining_seq;
			if ($str_len <= $max_str_len) {
				$keep_looping = 0;
				$aa_seq .= "\r\n" . $remaining_seq;
			} else {
				$aa_seq .= "\r\n" . substr($remaining_seq, 0, $max_str_len);
				$remaining_seq = substr($remaining_seq, $max_str_len) ;
			}
		}

		my $output_lines = $fasta_id . $aa_seq . "\r\n";

		open( TEMPFILE, ">$temp_file") or
			die "Cannot open $temp_file for writing\n";
		print TEMPFILE $output_lines;
		close TEMPFILE;

		# now run the tmhmm command, and capture the output in this program

		my $tmap_executable = '/home/emma/Emmas_files_on_linuxbox/software_to_install/emboss/EMBOSS-6.4.0/emboss/tmap';
		my $tmap_command = $tmap_executable . ' ' . $temp_file . ' -graph none -outfile ' . $temp_file_2;
		my $tmap_output = `$tmap_command`;

		# now process the tmap output

#		#=======================================
#		#
#		# Sequence: Consensus     from: 1   to: 417
#		# HitCount: 10
#		#=======================================
#
#		  Start     End TransMem Sequence
#		      8      36        1 NFWMFGLFFFFYFFIMGAYFPFFPIWLHD
#		     41      69        2 SKSDTGIIFAAISLFSLLFQPLFGLLSDK
#		     77     105        3 LWIITGMLVMFAPFFIFIFGPLLQYNILV
#		    144     172        4 RMFGCVGWALGASIVGIMFTINNQFVFWL
#		    219     237        5 QPKLWFLSLYVIGVSCTYD
#		    262     282        6 GYVTTMGELLNASIMFFAPLI
#		    288     308        7 GKNALLLAGTIMSVRIIGSSF
#		    314     334        8 EVVILKTLHMFEVPFLLVGCF
#		    344     372        9 RFSATIYLVCFCFFKQLAMIFMSVLAGNM
#		    375     403       10 SIGFQGAYLVLGLVALGFTLISVFTLSGP
#
#		#---------------------------------------
#		#---------------------------------------
#		#=======================================
#		#
#		# Sequence: 2cfp_A     from: 1   to: 417
#		# HitCount: 10
#		#=======================================
#
#		  Start     End TransMem Sequence
#		      8      36        1 NFWMFGLFFFFYFFIMGAYFPFFPIWLHD
#		     41      69        2 SKSDTGIIFAAISLFSLLFQPLFGLLSDK
#		     77     105        3 LWIITGMLVMFAPFFIFIFGPLLQYNILV
#		    144     172        4 RMFGCVGWALGASIVGIMFTINNQFVFWL
#		    219     237        5 QPKLWFLSLYVIGVSCTYD
#		    262     282        6 GYVTTMGELLNASIMFFAPLI
#		    288     308        7 GKNALLLAGTIMSVRIIGSSF
#		    314     334        8 EVVILKTLHMFEVPFLLVGCF
#		    344     372        9 RFSATIYLVCFCFFKQLAMIFMSVLAGNM
#		    375     403       10 SIGFQGAYLVLGLVALGFTLISVFTLSGP
#
#		#---------------------------------------

		open TMAP, $temp_file_2 or die $!;
		my @tmap_lines = <TMAP>;

		my $output_segments = '';

		my $start_reading_transmem = 0;
		my $start_of_transmem_text = '# Sequence: ' . $seqid;
		my $length_start_text = length($start_of_transmem_text);
		my $end_of_transmem_text = '#---------------------------------------';
		foreach my $tmap_line (@tmap_lines) {
			chomp($tmap_line);
			$tmap_line = trim($tmap_line);
			if ($tmap_line ne '') {
				if ($start_reading_transmem == -1) {
					my $do_nothing = 1;
				} elsif ($start_reading_transmem == 0) {
					if (length($tmap_line) >= $length_start_text) {
						if (substr($tmap_line,0,$length_start_text) eq $start_of_transmem_text) {
							$start_reading_transmem++;
						}
					}
				} elsif ($start_reading_transmem > 0) {
					if ($tmap_line eq $end_of_transmem_text) {
						$start_reading_transmem = -1;
					} else {
						if ($start_reading_transmem >= 4) {
							my $this_from = substr($tmap_line,0,7);
							my $this_to = substr($tmap_line,7,8);
							$this_from = trim($this_from);
							$this_to = trim($this_to);
							$output_segments .= $this_from . '-' . $this_to . ',';
						}
						$start_reading_transmem++;
					}
				}
			}
		}
		if ($output_segments eq '') {
			$output_segments = '-';
		}

		my $output_line = $seqid . ':::::' . $input_aaseq . ':::::' . $output_segments . ':::::';

		print OUTFILE "$output_line\n";
		close OUTFILE;

		open( OUTFILE, ">>$output_file") or
			die "Cannot open $output_file for appending\n";
	}
}
close OUTFILE;



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


