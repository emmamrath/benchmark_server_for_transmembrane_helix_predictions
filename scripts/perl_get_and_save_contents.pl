#!/usr/bin/perl   -w

# perl perl_get_and_save_contents.pl -outfile mempype_static_jquery.js"
# perl perl_get_and_save_contents.pl -outfile mempype_static_calendar.js"

# If the -rmvggg flag is present, then remove leading and trailing GGG, HHH and AAA before submitting to prediction tool.


use strict;
use Getopt::Long;
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
#my $url_for_get = 'http://mu2py.biocomp.unibo.it/mempype/static/jquery.js';
my $url_for_get = 'http://mu2py.biocomp.unibo.it/mempype/static/calendar.js';

my $max_num_retries_of_query1 = 10;
my $referer = 'Mozilla/5.0 (X11; Ubuntu; Linux i686; rv:10.0.2) Gecko/20100101 Firefox/10.0.2';

my $output_file;
GetOptions( "outfile=s" => \$output_file );
if (!defined $output_file) {
	die "Usage $0 -outfile OUTPUTfile\n";
}

# Open output file

open( OUTFILE, ">$output_file") or
	die "Cannot open $output_file for writing\n";

my $ua = LWP::UserAgent->new;
#my $newagent = 'test search program (' . $ua->agent . ')';
my $newagent = 'Mozilla/5.0 (X11; Ubuntu; Linux i686; rv:10.0.2) Gecko/20100101 Firefox/10.0.2';
$ua->agent($newagent);

my $request = HTTP::Request->new( GET => $url_for_get );

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

	$num_retries += 1;
}

if ($got_timeout == 0) {

	my $content = $response->content;

	my @lines = split(/\n/, $content);

	foreach my $line (@lines) {

		print OUTFILE "$line\n";
	}
}

close OUTFILE;



sub trim {
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}

