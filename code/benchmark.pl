#!/usr/bin/perl   -wT
use cPanelUserConfig;

# Copyright 2013 Emma Rath. All rights reserved.
# Soon this program will be released under an open source software license such as GNU General Public License or 
# Creative Commons license for Free Software Foundation's GNU General Public License at creativecommons.org

# INTRODUCTION :
# 
# This benchmark server allows you to compare membrane helix prediction and topology prediction methods that have made their predictions based on amino acid sequence.
# 
# To compare the various methods, choose the benchmark dataset including any restrictions to the dataset (or leave the choices as the default ones),
# and then press the RUN button further down below on this web page.
# 
# If you have your own membrane helix prediction method that you would like to benchmark against other methods,
# then click the 'display FASTA sequences' option before pressing the RUN button below.
# This will display the FASTA sequences for the benchmark dataset that you have chosen.
# After running those sequences on your method, upload your prediction results below before pressing the RUN button
# (this time with the default 'run the benchmark' option selected).
# 
# The membrane helices from 3D protein structures used to benchmark the various prediction methods have been taken from the RCSB Protein Data Bank (PDB).
# The membrane position predictions used on those proteins were taken from
# Orientations of Proteins in Membranes (OPM) database and PDBTM: Protein Data Bank of Transmembrane Proteins.

# If you use this server, please cite :
#
#	Emma M. Rath, Dominique Tessier, Hong Ching Lee, Tim Werner, Alexander A. Campbell, Noeris K. Salam, Lawrence K. Lee, W. Bret Church (2013)
#	'A benchmark server using high resolution protein structure data, and benchmark results for membrane helix predictions'
#	BMC Bioinformatics, 14, 111 doi: 10.1186/1471-2105-14-111.


# Reference :
#
#	Emma M. Rath, Dominique Tessier, Hong Ching Lee, Tim Werner, Alexander A. Campbell, Noeris K. Salam, Lawrence K. Lee, W. Bret Church (2013)
#	'A benchmark server using high resolution protein structure data, and benchmark results for membrane helix predictions'
#	BMC Bioinformatics, 14, 111 doi: 10.1186/1471-2105-14-111.
#
#	Kernytsky A, Rost B (2003)
#	'Static benchmarking of membrane helix predictions'
#	Nucleic Acids Research, 31, 3642-3644.
#
#	Chen CP, Kernytsky A, Rost B (2003)
#	'Transmembrane helix predictions revisited'
#	Protein Science, 11, 2774-2791.
#
#	Lomize MA, Lomize AL, Pogozheva ID, Mosberg HI (2006)
#	'OPM: Orientations of Proteins in Membranes database'
#	Bioinformatics,22, 623-625.
#
#	Kabsch W, Sander C (1983)
#	'Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features'
#	Biopolymers, 22, 2577-637.
#
#	Hobohm U, Scharf M, Schneider R, Sander C (1992)
#	'Selection of representative protein data sets'
#	Protein Science, 1, 409-417.
#
#	Sander C, Schneider R (1991)
#	'Database of Homology-Derived Protein Structures and the Structural Meaning of Sequence Alignment'
#	Proteins: Structure, Function, and Genetics, 9, 56-68.
#
#	Rice P, Longden I, Bleasby A (2000)
#	'EMBOSS: the European Molecular Biology open software suite'
#	Trends in Genetics, 16, 276-277.
#
#	Needleman SB, Wunsch CD (1970)
#	'A general method applicable to the search for similarities in the amino acid sequence of two proteins'
#	Journal of Molecular Biology, 48, 443-453.
#
#	Smith TF, Waterman MS (1981)
#	'Identification of Common Molecular Subsequences'
#	Journal of Molecular Biology, 147, 195-197.
#
#	Tusnády GE, Dosztányi Z, Simon I. (2004)
#	'Transmembrane proteins in the Protein Data Bank: identification and classification.'
#	Bioinformatics, 20, 2964-2972.
#
#	Tusnády GE, Dosztányi Z, Simon I. (2005)
#	'PDB_TM: selection and membrane localization of transmembrane proteins in the protein data bank.'
#	Nucleic Acids Research, 33(Database issue), D275-8.



use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use Math::Complex;
use Math::Trig;
use CGI;
use GD;
use GD::Text;
use GD::Graph::hbars;
use Data::Dumper;

#===================================================================================================
#	main
#===================================================================================================

my $q = new CGI;
my $global;
my $data;
my $output;
my $debug = '';

init();

get_program_mode();

if ($global->{'program_mode'} eq 'cgi-display-form') {

	display_html_form();

} elsif ($global->{'program_mode'} eq 'cgi-output') {

	process_input_produce_output_results();

} elsif ($global->{'program_mode'} eq 'cgi-output-graph') {

	display_html_graph();

}

#===================================================================================================
#	figure out whether this program is being called to display the HTML input form
#	or whether the form has been filled and now the results need to be displayed
#===================================================================================================

sub get_program_mode {

	my $params = $q->Vars;

	$global->{'program_mode'} = 'cgi-display-form';

	if (defined($params->{'graph_coords'})) {

		$global->{'program_mode'} = 'cgi-output-graph';

	} else {

		if (defined($ENV{'REQUEST_METHOD'})) {
			if ($ENV{'REQUEST_METHOD'} eq 'POST') {

				$global->{'program_mode'} = 'cgi-output';
			}
		}
	}
}

#===================================================================================================
#	the form has been filled and now the results need to be calculated and displayed
#===================================================================================================

sub process_input_produce_output_results {

	validate_and_build_results();

	if ($global->{'got_error'} == 1) {

		display_html_form_with_error();

	} else {

		if ($global->{'display_text_not_html'} == 1) {

			display_text_results();

		} elsif ($global->{'display_text_not_html'} == 2) {

			display_text_results();

		} else {

			if ($global->{'run'} eq 'benchmark') {

				display_html_benchmark_results();

			} elsif (($global->{'run'} eq 'display') or ($global->{'run'} eq 'display_date')) {

				display_html_aligned_sequences();

			} else {

				$global->{'got_error'} = 1;
				$global->{'err_msg'} .= "It seems that input specifying what to do was not received.<br>This is probably a system error.<br>\n";
				display_html_form_with_error();
			}
		}
	}
}

#===================================================================================================
#	validate_and_build_results
#===================================================================================================

sub validate_and_build_results {

	read_and_validate_input_parameters();

	if ($global->{'got_error'} == 0) {
		if ($global->{'manual_list_seqs'} eq '') {

			build_list_of_output_ids();

		}
		# if ($global->{'manual_list_seqs'} ne '') then we have already built the list of output ids when validating them.
	}

	validate_user_predictions_list();

	if ($global->{'got_error'} == 0) {

		if (($global->{'display_text_not_html'} == 1) or ($global->{'display_text_not_html'} == 2)) {
			if ( defined(%{$output->{'id'}->{'all'}}) ) { # if the list of valid sequences is not empty, process each sequence for output
				build_text_results();
			}

		} elsif ($global->{'run'} eq 'benchmark') {
			if ( defined(%{$output->{'id'}->{'all'}}) ) { # if the list of valid sequences is not empty, process each sequence for output
				build_benchmark_results();
			}

		} elsif (($global->{'run'} eq 'display') or ($global->{'run'} eq 'display_date')) {
			if ( defined(%{$output->{'id'}->{'all'}}) ) { # if the list of valid sequences is not empty, process each sequence for output
				build_view_results();
			}
		}
	}
}

#===================================================================================================
#	display_html_form
#===================================================================================================

sub display_html_form {

	print $q->header();

	my $html_file = 'benchmark.html';
	open HTMLFILE, $html_file or die $!;
	my @input_lines = <HTMLFILE>;

	foreach my $input_line (@input_lines) {
		chomp $input_line;
		if ($input_line ne '') {

			my $output_line = $input_line . "\n";
			print $output_line;
		}
	}

	if ($debug ne '') {
		print "<br>\n$debug<br>\n";
	}
}

#===================================================================================================
#	display_html_form_with_error
#===================================================================================================

sub display_html_form_with_error {

	print $q->header();

	my $html_file = 'benchmark.html';
	open HTMLFILE, $html_file or die $!;
	my @input_lines = <HTMLFILE>;
	foreach my $input_line (@input_lines) {
		chomp $input_line;
		if ($input_line ne '') {
			my $output_line = $input_line . "\n";
			print $output_line;

			if ($input_line eq "<a name='msg'>") {
				$output_line = "<span class='red_bold'>" . $global->{'err_msg'} . "</span><br><br>\n\n";
				print $output_line;
			}
		}
	}

	if ($debug ne '') {
		print "<br>\n$debug<br>\n";
	}
}

#===================================================================================================
#	display_text_results
#===================================================================================================

sub display_text_results {

	# user requested to display just the protein-chain ids, or with FASTA sequence, or with protein name
	# if (($run eq 'fasta') or ($run eq 'list') or ($run eq 'names')) {

	if ($global->{'display_text_not_html'} == 2) {

		print $q->header("Content-type: text/plain"); # this instructs the user browser to save a download file consisting of the output from this program

	} else { # ($global->{'display_text_not_html'} == 1)

		print $q->header( -type => 'text/plain' ); # this instructs the user browser to display the output from this program on the screen as text, without needing html formatting
	}

	my $got_output_lines = 0;

	if ( defined(%{$output->{'id'}->{'all'}}) ) { # if the list of valid sequences is not empty, process each sequence for output

		if ($global->{'run'} eq 'list') {

			my %output_id = %{$output->{'id'}->{'all'}};

			foreach my $this_id ( sort keys %output_id ) {
				my $output_line = "$this_id\n";
				print $output_line;
				$got_output_lines = 1;
			}
		}

		if ($global->{'run'} eq 'names') {

			my %output_id = %{$output->{'id'}->{'all'}};
			my %output_name = %{$output->{'name'}};
			my %output_release_dates = %{$output->{'release_dates'}};
			my %output_pdb_methods = %{$output->{'pdb_methods'}};
			my %output_pdb_resolutions = %{$output->{'pdb_resolutions'}};
			my %output_mpstruc_category = %{$output->{'mpstruc_category'}};

			print "ID     RELEASE-DATE   EXPERIMENTAL-METHOD  RESOLUTION    NAME                                                                                                                                                               FAMILY\n";
			print "==     ============   ===================  ==========    ====                                                                                                                                                               ======\n";
			foreach my $this_id ( sort keys %output_id ) {
				my $this_name = $output_name{$this_id};
				$this_name = sprintf("%-160s", $this_name);
				my $this_release_date = $output_release_dates{$this_id};
				my $this_pdb_method = $output_pdb_methods{$this_id};
				$this_pdb_method = sprintf("%-24s", $this_pdb_method);
				my $this_pdb_resolution = '     ';
				if (exists($output_pdb_resolutions{$this_id})) {
					$this_pdb_resolution = $output_pdb_resolutions{$this_id};
					$this_pdb_resolution = sprintf("%-5s", $this_pdb_resolution);
				}
				my $this_mpstruc_category = $output_mpstruc_category{$this_id};
				my $output_line = "$this_id   $this_release_date   $this_pdb_method   $this_pdb_resolution   $this_name   $this_mpstruc_category\n";
				print $output_line;
				$got_output_lines = 1;
			}
		}

		if (($global->{'run'} eq 'fasta') or ($global->{'run'} eq 'fastaG') or ($global->{'run'} eq 'dload') or ($global->{'run'} eq 'dloadG')) {

			my %output_id = %{$output->{'id'}->{'all'}};
			my %output_aaseq = %{$output->{'aaseq'}};
			my $line_length = 60;
			my $line_length_minus_1 = $line_length - 1;

			foreach my $this_id ( sort keys %output_id ) {
				my $this_aaseq = $output_aaseq{$this_id};
				if (($global->{'run'} eq 'fastaG') or ($global->{'run'} eq 'dloadG')) {
					$this_aaseq =~ s/\_/G/g;
				}
				my $output_line = ">$this_id\n";
				for ( my $i = 0; $i < length($this_aaseq); $i++ ) {
					my $this_char = substr( $this_aaseq, $i, 1 );
					$output_line .= $this_char;
					if ( ($i % $line_length) == $line_length_minus_1 ) {
						$output_line .= "\n";
					}
				}

				if ( (length($this_aaseq) % $line_length) != 0 ) {
					$output_line .= "\n";
				}
				print $output_line;
				$got_output_lines = 1;
			}
		}

		if (($global->{'run'} eq 'list') or ($global->{'run'} eq 'names') or ($global->{'run'} eq 'fasta') or ($global->{'run'} eq 'fastaG') or ($global->{'run'} eq 'dload') or ($global->{'run'} eq 'dloadG')) {
			if ($got_output_lines == 0) {
				my $output_line = "No sequences are in requested list of sequences.\n";
				print $output_line;
			}
		}

	} else {
		my $output_line = "No valid sequences are in requested list of sequences.\n";
		print $output_line;
	}

	if ($global->{'run'} eq 'count') {

		my $count = 0;

		if ( defined(%{$output->{'id'}->{'all'}}) ) { # if the list of valid sequences is not empty, process each sequence for output

			my %output_id = %{$output->{'id'}->{'all'}};

			foreach my $this_id ( keys %output_id ) {
				$count++;
			}
		}

		my $output_line = "Number of sequences = $count\n";
		print $output_line;
			$got_output_lines = 1;
	}

	if ($debug ne '') {
		print "\n$debug\n";
	}
}

#===================================================================================================
#	display_html_graph
#===================================================================================================

sub display_html_graph {

	my $params = $q->Vars;
	my $graph_coords = '';
	if (defined($params->{'graph_coords'})) {
		$graph_coords = trim($params->{'graph_coords'});
	}
	my $graph_title = '';
	if (defined($params->{'graph_title'})) {
		$graph_title = trim($params->{'graph_title'});
	}
	my $graph_legend = '';
	if (defined($params->{'graph_legend'})) {
		$graph_legend = trim($params->{'graph_legend'});
	}

	#my @graph_data_2 = (
	#["MemBrain","VALPRED","HMMTOP2","TMLOOP","DASTMfilt","SVMtop","Phobius", "Philius", "PRED-TMR"],
	#[ 5, 12, 24, 33, 19, 8, 6, 15, 21],
	#[ 1, 2, 5, 6, 3, 1.5, 1, 3, 4],
	#[ 91, 82, 85, 96, 93, 915, 81, 83, 74],
	#); 
	my $num_methods = 0;
	my @graph_data;
	my @bits = split( /:::::/, $graph_coords );
	my $i = 0;
	foreach my $this_data_set (@bits) {
		if ($i == 0) {
			my @method_names;
			my @bits2 = split( /;;;;;/, $this_data_set );
			foreach my $this_method (@bits2) {
				push( @method_names, $this_method );
				$num_methods++;
			}
			push( @graph_data, \@method_names );
		} else {
			my @method_scores;
			my @bits2 = split( /;;;;;/, $this_data_set );
			foreach my $this_score (@bits2) {
				push( @method_scores, $this_score );
			}
			push( @graph_data, \@method_scores );
		}
		$i++;
	}

	# $graph_object->set_legend('Qhtm obs%', 'Qhtm prd%', 'Qok%');
	my @data_legend;
	$graph_legend =~ s/&#37;/%/g;
	my @bits3 = split( /:::::/, $graph_legend );
	foreach my $this_criteria (@bits3) {
		push( @data_legend, $this_criteria );
	}
	push( @data_legend, 'two' );
	push( @data_legend, 'three' );

	my $graph_width = 600;
	my $graph_height = 10 * $num_methods;
	if ($graph_height < 800) {
		$graph_height = 800;
	}

	my $graph_object = new GD::Graph::hbars( $graph_width, $graph_height );

	$graph_object->set(
	x_label => 'Methods',
	y_label => 'Benchmark Criteria',
	title => $graph_title,
	long_ticks => 1,
	y_max_value => 100,
	y_tick_number => 10,
	y_label_skip => 2,
	bar_spacing => 10,
	bar_width => 20,
	#shadow_depth => 4,
	accent_threshold => 200,
	transparent => 0,
	);

	$graph_object->set_title_font(['arial', 'verdana', gdGiantFont], 12);
	$graph_object->set_x_label_font(['arial', 'verdana', gdGiantFont], 12);
	$graph_object->set_y_label_font(['arial', 'verdana', gdGiantFont], 12);
	$graph_object->set_x_axis_font(['arial', 'verdana', gdGiantFont], 12);
	$graph_object->set_y_axis_font(['arial', 'verdana', gdGiantFont], 12);
	$graph_object->set_legend_font(['arial', 'verdana', gdGiantFont], 12);
	$graph_object->set_legend( @data_legend );
	$graph_object->plot(\@graph_data);

	# Apparently we need to do this on Windows or the image will be garbled, and it doesn't hurt on Unix/Linux/etc.  
	binmode STDOUT;

	my $ext = $graph_object->export_format;
	print $q->header("Content-type: image/gif");
	print $graph_object->gd->$ext(); 
}

#===================================================================================================
#	display_html_benchmark_results
#===================================================================================================

sub display_html_benchmark_results {

	print $q->header();

	#===================================================================================================
	#	generic top of the displayed html page
	#===================================================================================================

	my $html_file = 'benchmark.html';
	open HTMLFILE, $html_file or die $!;
	my @input_lines = <HTMLFILE>;

	my $stop_writing_html_from_file = 0;
	foreach my $input_line (@input_lines) {

		if ($stop_writing_html_from_file == 0) {

			chomp $input_line;
			if ($input_line ne '') {

				my $output_line = $input_line . "\n";
				print $output_line;

				if ($input_line eq "<a name='msg'>") {
					$stop_writing_html_from_file = 1;
				}
			}
		}
	}

	if ($global->{'info_msg'} ne '') {
		my $output_line = "<span class='red_bold'>" . $global->{'info_msg'} . "</span><br>\n\n";
		print $output_line;
	}

	#===================================================================================================
	#	the displayed benchmark results - summary
	#===================================================================================================

	print "<br>\n";

	print "<table cellspacing=0 cellpadding=10 border=1 bgcolor='#eeeeee' width=800><tr><td valign='top' align='left'>\n";
	print "<table cellspacing=0 cellpadding=0 border=0>\n";
	print "<tr><td valign='top' align='left'><span class='h1'>PARAMETERS AND DATA VOLUME USED TO CALCULATE BENCHMARK RESULTS :</span><br><br>\n";

	print "Number of unique benchmark sequences : " . $output->{'num_reference_seqs'} . "<br>\n";
	print "Number of benchmark membrane helices in those sequences : " . $output->{'num_reference_membrane_helices'} . "<br><br>\n";
	if ($global->{'use_the_user_predictions'} == 1) {
		print "Number of sequences in uploaded predictions data : " . $output->{'num_user_seqs'} . "<br>\n";
		print "Number of membrane helices in those uploaded sequences : " . $output->{'num_user_membrane_helices'} . "<br><br>\n";
	}
	if ($global->{'sort_by_name'} ne '') {
		print "Results sorted by column : " . $global->{'sort_by_name'} . "<br><br>\n";
	}
	print "Minimum length of observed and predicted helices (smaller helices are ignored) : " . $global->{'minimum_helix_length'} . "<br>\n";
	print "Minimum overlap of observed helix residues for a helix prediction to score in the<br>&nbsp;&nbsp;&nbsp;per-sequence-accuracy, per-segment-accuracy and average-helix-end-difference scores : " . $global->{'minimum_overlap'} . "<br>\n";
	print "Maximum distance in residues to include in the membrane helix approximate-ends scores : " . $global->{'max_dist_approx_ends'} . "<br>\n";
	print "<br>\n";

	if ($global->{'user_predictions_missing_defined_sequences'} == 1) {
		print "<p><font color='red'>The uploaded predictions data did not contain all the sequences defined in the chosen benchmark dataset.<br>This will adversely affect the benchmark statistics for the YOU method.</b><br><br></font></p>\n",
	}
	if ($global->{'user_predictions_have_extra_sequences'} == 1) {
		print "<p><font color='red'>The uploaded predictions data contained extra sequences than those defined in the chosen benchmark dataset.<br>The extra sequences were ignored and are not included in the benchmark statistics.</b><br><br></font></p>\n",
	}
	if ($global->{'user_file_and_textbox_predictions_present'} == 1) {
		print "<p><font color='red'>Two sets of predictions were entered, in an upload file and in the text box. The upload file was used and the text box entries were ignored.</b><br><br></font></p>\n",
	}
	if ($global->{'minimum_overlap_not_valid'} == 1) {
		print "<p><font color='red'>The entered 'Minimum overlap' parameter was not valid, so " . $data->{'minimum_overlap'} . " was used instead.</b><br><br></font></p>\n",
	}
	if ($global->{'minimum_helix_length_not_valid'} == 1) {
		print "<p><font color='red'>The entered 'Minimum length of observed and predicted helices' parameter was not valid, so " . $data->{'minimum_helix_length'} . " was used instead.</b><br><br></font></p>\n",
	}
	if ($global->{'max_dist_approx_ends_not_valid'} == 1) {
		print "<p><font color='red'>The entered 'Maximum distance in residues to include in the membrane-helix approximate-ends scores' parameter was not valid, so " . $data->{'max_dist_approx_ends'} . " was used instead.</b><br><br></font></p>\n",
	}
	print "</td></tr></table></td></tr></table><br>\n";

	#===================================================================================================
	#	display the relevent legend at the top of the page before showing the results
	#===================================================================================================

	my $start_writing_html_from_file = 0;
	foreach my $input_line (@input_lines) {

		chomp $input_line;
		if ($input_line ne '') {

			if ($input_line eq "<div id='div_for_key_to_benchmark_results_show'>") {
				$start_writing_html_from_file = 1;
			}
			elsif ($input_line eq "<a name='key_view'>") {
				$start_writing_html_from_file = 0;
			}
			elsif ($input_line eq "<div id='div_for_key_to_display_aligned_show'>") {
				$start_writing_html_from_file = 0;
			}

			if ($start_writing_html_from_file == 1) {
				my $output_line = $input_line . "\n";
				print $output_line;
			}
		}
	}
	my $output_line = "<br><br>\n\n";
	print $output_line;

	#===================================================================================================
	#	the displayed benchmark results - details
	#===================================================================================================

	if ( defined(@{$output->{'qok'}}) ) { # if the list of valid sequences is not empty, process each sequence for output

		my @output_method = @{$output->{'prediction_method'}};
		my @output_qok = @{$output->{'qok'}};
		my @output_qhtm_obs = @{$output->{'qhtm_obs'}};
		my @output_qhtm_prd = @{$output->{'qhtm_prd'}};
		my @output_avhe_diff = @{$output->{'avhe_diff'}};
		my @output_qhe_obs = @{$output->{'qhe_obs'}};
		my @output_qhe_prd = @{$output->{'qhe_prd'}};
		my @output_qnhe_obs = @{$output->{'qnhe_obs'}};
		my @output_qnhe_prd = @{$output->{'qnhe_prd'}};
		my @output_q2 = @{$output->{'q2'}};
		my @output_htm_mcc = @{$output->{'htm_mcc'}};
		my @output_q2t_obs = @{$output->{'q2t_obs'}};
		my @output_q2t_prd = @{$output->{'q2t_prd'}};
		my @output_q2n_obs = @{$output->{'q2n_obs'}};
		my @output_q2n_prd = @{$output->{'q2n_prd'}};
		my @output_qok3 = @{$output->{'qok3'}};
		my @output_nterm = @{$output->{'nterm'}};
		my @output_q2_ioseg = @{$output->{'q2_ioseg'}};
		my @output_qiom_obs = @{$output->{'qiom_obs'}};
		my @output_qiom_prd = @{$output->{'qiom_prd'}};
		my @output_qio_obs = @{$output->{'qio_obs'}};
		my @output_qio_prd = @{$output->{'qio_prd'}};
		my @output_q3 = @{$output->{'q3'}};
		my @output_q2_iores = @{$output->{'q2_iores'}};
		my @output_io_mcc = @{$output->{'io_mcc'}};

		my $graph_coords = '';
		my $graph_legend = '';
		my $graph_title = $global->{'graph_title'};
		for ( my $i = 0; $i < @output_method; $i++ ) {
			my $method = $output_method[$i];
			my $method_name = $data->{'method_short_name'}->{$method};
			$graph_coords .= $method_name . ';;;;;';
		}
		$graph_coords .= ':::::';
		my @list_of_benchmark_criteria = @{$data->{'list_of_benchmark_criteria'}};
		my @list_of_benchmark_criteria_name = @{$data->{'list_of_benchmark_criteria_name'}};
		if ($global->{'graph_by'} ne '') {
			my $this_graph_by = $global->{'graph_by'};
			my @output_criteria;
			if ($this_graph_by eq 'qok') {
				@output_criteria = @output_qok;
			} elsif ($this_graph_by eq 'qhtm_obs') {
				@output_criteria = @output_qhtm_obs;
			} elsif ($this_graph_by eq 'qhtm_prd') {
				@output_criteria = @output_qhtm_prd;
			} elsif ($this_graph_by eq 'avhe_diff') {
				@output_criteria = @output_avhe_diff;
			} elsif ($this_graph_by eq 'qhe_obs') {
				@output_criteria = @output_qhe_obs;
			} elsif ($this_graph_by eq 'qhe_prd') {
				@output_criteria = @output_qhe_prd;
			} elsif ($this_graph_by eq 'qnhe_obs') {
				@output_criteria = @output_qnhe_obs;
			} elsif ($this_graph_by eq 'qnhe_prd') {
				@output_criteria = @output_qnhe_prd;
			} elsif ($this_graph_by eq 'q2') {
				@output_criteria = @output_q2;
			} elsif ($this_graph_by eq 'htm_mcc') {
				@output_criteria = @output_htm_mcc;
			} elsif ($this_graph_by eq 'q2t_obs') {
				@output_criteria = @output_q2t_obs;
			} elsif ($this_graph_by eq 'q2t_prd') {
				@output_criteria = @output_q2t_prd;
			} elsif ($this_graph_by eq 'q2n_obs') {
				@output_criteria = @output_q2n_obs;
			} elsif ($this_graph_by eq 'q2n_prd') {
				@output_criteria = @output_q2n_prd;
			} elsif ($this_graph_by eq 'qok3') {
				@output_criteria = @output_qok3;
			} elsif ($this_graph_by eq 'nterm') {
				@output_criteria = @output_nterm;
			} elsif ($this_graph_by eq 'q2_ioseg') {
				@output_criteria = @output_q2_ioseg;
			} elsif ($this_graph_by eq 'qiom_obs') {
				@output_criteria = @output_qiom_obs;
			} elsif ($this_graph_by eq 'qiom_prd') {
				@output_criteria = @output_qiom_prd;
			} elsif ($this_graph_by eq 'qio_obs') {
				@output_criteria = @output_qio_obs;
			} elsif ($this_graph_by eq 'qio_prd') {
				@output_criteria = @output_qio_prd;
			} elsif ($this_graph_by eq 'q3') {
				@output_criteria = @output_q3;
			} elsif ($this_graph_by eq 'q2_iores') {
				@output_criteria = @output_q2_iores;
			} elsif ($this_graph_by eq 'io_mcc') {
				@output_criteria = @output_io_mcc;
			}
			for ( my $i = 0; $i < @output_method; $i++ ) {
				my $method = $output_method[$i];
				$graph_coords .= $output_criteria[$i] . ';;;;;';
			}
			$graph_coords .= ':::::';
			my $i = 0;
			foreach my $this_criteria (@list_of_benchmark_criteria) {
				if ($this_criteria eq $this_graph_by) {
					$graph_legend .= $list_of_benchmark_criteria_name[$i] . ':::::';
				}
				$i++;
			}
		}
		if ($global->{'graph_by_2'} ne '') {
			my $this_graph_by = $global->{'graph_by_2'};
			my @output_criteria;
			if ($this_graph_by eq 'qok') {
				@output_criteria = @output_qok;
			} elsif ($this_graph_by eq 'qhtm_obs') {
				@output_criteria = @output_qhtm_obs;
			} elsif ($this_graph_by eq 'qhtm_prd') {
				@output_criteria = @output_qhtm_prd;
			} elsif ($this_graph_by eq 'avhe_diff') {
				@output_criteria = @output_avhe_diff;
			} elsif ($this_graph_by eq 'qhe_obs') {
				@output_criteria = @output_qhe_obs;
			} elsif ($this_graph_by eq 'qhe_prd') {
				@output_criteria = @output_qhe_prd;
			} elsif ($this_graph_by eq 'qnhe_obs') {
				@output_criteria = @output_qnhe_obs;
			} elsif ($this_graph_by eq 'qnhe_prd') {
				@output_criteria = @output_qnhe_prd;
			} elsif ($this_graph_by eq 'q2') {
				@output_criteria = @output_q2;
			} elsif ($this_graph_by eq 'htm_mcc') {
				@output_criteria = @output_htm_mcc;
			} elsif ($this_graph_by eq 'q2t_obs') {
				@output_criteria = @output_q2t_obs;
			} elsif ($this_graph_by eq 'q2t_prd') {
				@output_criteria = @output_q2t_prd;
			} elsif ($this_graph_by eq 'q2n_obs') {
				@output_criteria = @output_q2n_obs;
			} elsif ($this_graph_by eq 'q2n_prd') {
				@output_criteria = @output_q2n_prd;
			} elsif ($this_graph_by eq 'qok3') {
				@output_criteria = @output_qok3;
			} elsif ($this_graph_by eq 'nterm') {
				@output_criteria = @output_nterm;
			} elsif ($this_graph_by eq 'q2_ioseg') {
				@output_criteria = @output_q2_ioseg;
			} elsif ($this_graph_by eq 'qiom_obs') {
				@output_criteria = @output_qiom_obs;
			} elsif ($this_graph_by eq 'qiom_prd') {
				@output_criteria = @output_qiom_prd;
			} elsif ($this_graph_by eq 'qio_obs') {
				@output_criteria = @output_qio_obs;
			} elsif ($this_graph_by eq 'qio_prd') {
				@output_criteria = @output_qio_prd;
			} elsif ($this_graph_by eq 'q3') {
				@output_criteria = @output_q3;
			} elsif ($this_graph_by eq 'q2_iores') {
				@output_criteria = @output_q2_iores;
			} elsif ($this_graph_by eq 'io_mcc') {
				@output_criteria = @output_io_mcc;
			}
			for ( my $i = 0; $i < @output_method; $i++ ) {
				my $method = $output_method[$i];
				$graph_coords .= $output_criteria[$i] . ';;;;;';
			}
			$graph_coords .= ':::::';
			my $i = 0;
			foreach my $this_criteria (@list_of_benchmark_criteria) {
				if ($this_criteria eq $this_graph_by) {
					$graph_legend .= $list_of_benchmark_criteria_name[$i] . ':::::';
				}
				$i++;
			}
		}
		if ($global->{'graph_by_3'} ne '') {
			my $this_graph_by = $global->{'graph_by_3'};
			my @output_criteria;
			if ($this_graph_by eq 'qok') {
				@output_criteria = @output_qok;
			} elsif ($this_graph_by eq 'qhtm_obs') {
				@output_criteria = @output_qhtm_obs;
			} elsif ($this_graph_by eq 'qhtm_prd') {
				@output_criteria = @output_qhtm_prd;
			} elsif ($this_graph_by eq 'avhe_diff') {
				@output_criteria = @output_avhe_diff;
			} elsif ($this_graph_by eq 'qhe_obs') {
				@output_criteria = @output_qhe_obs;
			} elsif ($this_graph_by eq 'qhe_prd') {
				@output_criteria = @output_qhe_prd;
			} elsif ($this_graph_by eq 'qnhe_obs') {
				@output_criteria = @output_qnhe_obs;
			} elsif ($this_graph_by eq 'qnhe_prd') {
				@output_criteria = @output_qnhe_prd;
			} elsif ($this_graph_by eq 'q2') {
				@output_criteria = @output_q2;
			} elsif ($this_graph_by eq 'htm_mcc') {
				@output_criteria = @output_htm_mcc;
			} elsif ($this_graph_by eq 'q2t_obs') {
				@output_criteria = @output_q2t_obs;
			} elsif ($this_graph_by eq 'q2t_prd') {
				@output_criteria = @output_q2t_prd;
			} elsif ($this_graph_by eq 'q2n_obs') {
				@output_criteria = @output_q2n_obs;
			} elsif ($this_graph_by eq 'q2n_prd') {
				@output_criteria = @output_q2n_prd;
			} elsif ($this_graph_by eq 'qok3') {
				@output_criteria = @output_qok3;
			} elsif ($this_graph_by eq 'nterm') {
				@output_criteria = @output_nterm;
			} elsif ($this_graph_by eq 'q2_ioseg') {
				@output_criteria = @output_q2_ioseg;
			} elsif ($this_graph_by eq 'qiom_obs') {
				@output_criteria = @output_qiom_obs;
			} elsif ($this_graph_by eq 'qiom_prd') {
				@output_criteria = @output_qiom_prd;
			} elsif ($this_graph_by eq 'qio_obs') {
				@output_criteria = @output_qio_obs;
			} elsif ($this_graph_by eq 'qio_prd') {
				@output_criteria = @output_qio_prd;
			} elsif ($this_graph_by eq 'q3') {
				@output_criteria = @output_q3;
			} elsif ($this_graph_by eq 'q2_iores') {
				@output_criteria = @output_q2_iores;
			} elsif ($this_graph_by eq 'io_mcc') {
				@output_criteria = @output_io_mcc;
			}
			for ( my $i = 0; $i < @output_method; $i++ ) {
				my $method = $output_method[$i];
				$graph_coords .= $output_criteria[$i] . ';;;;;';
			}
			$graph_coords .= ':::::';
			my $i = 0;
			foreach my $this_criteria (@list_of_benchmark_criteria) {
				if ($this_criteria eq $this_graph_by) {
					$graph_legend .= $list_of_benchmark_criteria_name[$i] . ':::::';
				}
				$i++;
			}
		}

		my $bg_colour = '#eeeeee';
		#my $panel_colour = '#dddddd';
		my $panel_colour = '#eeeeee';
		my $topography_colour = '#c8c8ed';
		my $topology_colour = '#c8edd1';

		print "<table cellspacing=0 cellpadding=0 border=1>\n<tr>\n<td>\n";

		print "<table cellspacing=0 cellpadding=10 border=0>\n";

		print "<tr></td><td valign='top' align='left' colspan=29 bgcolor='$bg_colour'><span class='h1'>BENCHMARK RESULTS :</span></td></tr>\n";


		print "<tr>\n",

		"			<td valign='bottom' align='left' colspan=1 bgcolor='$panel_colour'><span class='bold'></span></td>\n",
		"			<td width=1 bgcolor='$bg_colour'></td>\n",
		"			<td valign='bottom' align='center' colspan=15 bgcolor='$topography_colour'><span class='bold'>...Membrane helices versus not membrane helix...</span></td>\n",
		"			<td width=1 bgcolor='$bg_colour'></td>\n",
		"			<td valign='bottom' align='center' colspan=10 bgcolor='$topology_colour'><span class='bold'>...Inside/outside topology...</span></td>\n",
		"			<td width=1 bgcolor='$bg_colour'></td>\n",
		"</tr>\n";

		print "<tr>\n",

		"			<td valign='bottom' align='left' colspan=1 bgcolor='$panel_colour'><span class='bold'>Method</span></td>\n",
		"			<td width=1 bgcolor='$bg_colour'></td>\n",

		"			<td valign='bottom' align='right' colspan=1 bgcolor='$topography_colour'><span class='bold'>per sequence</span></td>\n",
		"			<td width=1 bgcolor='$topography_colour'></td>\n",
		"			<td valign='bottom' align='left' colspan=2 bgcolor='$topography_colour'><span class='bold'>per segment</span></td>\n",
		"			<td width=1 bgcolor='$topography_colour'></td>\n",
		"			<td valign='bottom' align='left' colspan=3 bgcolor='$topography_colour'><span class='bold'>per helix end</span></td>\n",
		"			<td width=1 bgcolor='$topography_colour'></td>\n",
		"			<td valign='bottom' align='left' colspan=6 bgcolor='$topography_colour'><span class='bold'>per residue</span></td>\n",
		"			<td width=1 bgcolor='$bg_colour'></td>\n",

		"			<td valign='bottom' align='right' colspan=2 bgcolor='$topology_colour'><span class='bold'>per sequence</span></td>\n",
		"			<td width=1 bgcolor='$topology_colour'></td>\n",
		"			<td valign='bottom' align='left' colspan=3 bgcolor='$topology_colour'><span class='bold'>per segment</span></td>\n",
		"			<td width=1 bgcolor='$topology_colour'></td>\n",
		"			<td valign='bottom' align='left' colspan=3 bgcolor='$topology_colour'><span class='bold'>per residue</span></td>\n",
		"			<td width=1 bgcolor='$bg_colour'></td>\n",
		"</tr>\n";

		print "<tr>\n";
		print "			<td valign='bottom' align='left' bgcolor='$panel_colour'><b></span></td>\n";
		print "			<td width=1 bgcolor='$bg_colour'></td>\n";

		print "			<td valign='bottom' align='right' bgcolor='$topography_colour'><span class='bold'>Qok<br>&#37;</span></td>\n";
		print "			<td bgcolor='$topography_colour'></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topography_colour'><span class='bold'><span class='smaller_bold'>(sensitivity)</span><br>Qhtm<br>&#37;obs</span></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topography_colour'><span class='bold'><span class='smaller_bold'>(specificity)</span><br>Qhtm<br>&#37;prd</span></td>\n";
		print "			<td bgcolor='$topography_colour'></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topography_colour'><span class='bold'>AvHb diff</span><br><span class='smaller_bold'>in residues</span></span></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topography_colour'><span class='bold'>QHb<br>&#37;obs</span></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topography_colour'><span class='bold'>Gauss<br>QHb<br>&#37;obs</span></td>\n";
		print "			<td bgcolor='$topography_colour'></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topography_colour'><span class='bold'>Q2<br>&#37;</span></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topography_colour'><span class='bold'>htm<br>MCC<br><span class='smaller_bold'>(-1..0..1)</span></span></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topography_colour'><span class='bold'>Q2T<br>&#37;obs</span></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topography_colour'><span class='bold'>Q2T<br>&#37;prd</span></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topography_colour'><span class='bold'>Q2N<br>&#37;obs</span></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topography_colour'><span class='bold'>Q2N<br>&#37;prd</span></td>\n";
		print "			<td bgcolor='$bg_colour'></td>\n";

		print "			<td valign='bottom' align='right' bgcolor='$topology_colour'><span class='bold'>Qok3<br>&#37;</span></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topology_colour'><span class='bold'>Nterm<br>&#37;</span></td>\n";
		print "			<td bgcolor='$topology_colour'></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topology_colour'><span class='bold'>ioSeg<br>Q2<br>&#37;</span></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topology_colour'><span class='bold'>Qiom<br>&#37;obs</span></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topology_colour'><span class='bold'>Qio<br>&#37;obs</span></td>\n";
		print "			<td bgcolor='$topology_colour'></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topology_colour'><span class='bold'>Q3<br>&#37;</span></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topology_colour'><span class='bold'>ioRes<br>Q2<br>&#37;</span></td>\n";
		print "			<td valign='bottom' align='right' bgcolor='$topology_colour'><span class='bold'>io<br>MCC<br><span class='smaller_bold'>(-1..0..1)</span></span></td>\n";
		print "			<td bgcolor='$bg_colour'></td>\n";
		print "</tr>\n";

		for ( my $i = 0; $i < @output_method; $i++ ) {

			my $method = $output_method[$i];
			my $method_name = $data->{'method_long_name'}->{$method};

			my $qok = $output_qok[$i];
			my $qhtm_obs = $output_qhtm_obs[$i];
			my $qhtm_prd = $output_qhtm_prd[$i];
			my $avhe_diff = $output_avhe_diff[$i];
			my $qhe_obs = $output_qhe_obs[$i];
			my $qhe_prd = $output_qhe_prd[$i];
			my $qnhe_obs = $output_qnhe_obs[$i];
			my $qnhe_prd = $output_qnhe_prd[$i];
			my $q2 = $output_q2[$i];
			my $htm_mcc = $output_htm_mcc[$i];
			my $q2t_obs = $output_q2t_obs[$i];
			my $q2t_prd = $output_q2t_prd[$i];
			my $q2n_obs = $output_q2n_obs[$i];
			my $q2n_prd = $output_q2n_prd[$i];

			my $qok3 = $output_qok3[$i];
			my $nterm = $output_nterm[$i];
			my $q2_ioseg = $output_q2_ioseg[$i];
			my $qiom_obs = $output_qiom_obs[$i];
			my $qiom_prd = $output_qiom_prd[$i];
			my $qio_obs = $output_qio_obs[$i];
			my $qio_prd = $output_qio_prd[$i];
			my $q3 = $output_q3[$i];
			my $q2_iores = $output_q2_iores[$i];
			my $io_mcc = $output_io_mcc[$i];

			my $b_flag = "";
			my $b_flag2 = "";

			print "		<tr>\n";
			print "			<td valign='top' align='left' bgcolor='$bg_colour'><span class='smaller'>$b_flag$method_name$b_flag2</span></td>\n";
			print "			<td bgcolor='$bg_colour'></td>\n";

			print "			<td valign='top' align='right'>$b_flag$qok$b_flag2</td>\n";
			print "			<td bgcolor='$bg_colour'></td>\n";

			print "			<td valign='top' align='right'>$b_flag$qhtm_obs$b_flag2</td>\n";
			print "			<td valign='top' align='right'>$b_flag$qhtm_prd$b_flag2</td>\n";
			print "			<td bgcolor='$bg_colour'></td>\n";
			print "			<td valign='top' align='right'>$b_flag$avhe_diff$b_flag2</td>\n";
			print "			<td valign='top' align='right'>$b_flag$qhe_obs$b_flag2</td>\n";
			print "			<td valign='top' align='right'>$b_flag$qnhe_obs$b_flag2</td>\n";
			print "			<td bgcolor='$bg_colour'></td>\n";
			print "			<td valign='top' align='right'>$b_flag$q2$b_flag2</td>\n";
			print "			<td valign='top' align='right'>$b_flag$htm_mcc$b_flag2</td>\n";
			print "			<td valign='top' align='right'>$b_flag$q2t_obs$b_flag2</td>\n";
			print "			<td valign='top' align='right'>$b_flag$q2t_prd$b_flag2</td>\n";
			print "			<td valign='top' align='right'>$b_flag$q2n_obs$b_flag2</td>\n";
			print "			<td valign='top' align='right'>$b_flag$q2n_prd$b_flag2</td>\n";
			print "			<td bgcolor='$bg_colour'></td>\n";

			print "			<td valign='top' align='right'>$b_flag$qok3$b_flag2</td>\n";
			print "			<td valign='top' align='right'>$b_flag$nterm$b_flag2</td>\n";
			print "			<td bgcolor='$bg_colour'></td>\n";
			print "			<td valign='top' align='right'>$b_flag$q2_ioseg$b_flag2</td>\n";
			print "			<td valign='top' align='right'>$b_flag$qiom_obs$b_flag2</td>\n";
			print "			<td valign='top' align='right'>$b_flag$qio_obs$b_flag2</td>\n";
			print "			<td bgcolor='$bg_colour'></td>\n";
			print "			<td valign='top' align='right'>$b_flag$q3$b_flag2</td>\n";
			print "			<td valign='top' align='right'>$b_flag$q2_iores$b_flag2</td>\n";
			print "			<td valign='top' align='right'>$b_flag$io_mcc$b_flag2</td>\n";
			print "			<td bgcolor='$bg_colour'></td>\n";
			print "		</tr>\n";
		}

		print "<tr><td height=10 colspan=29 bgcolor='$bg_colour'></td></tr>\n";
		if ($global->{'kingdom'}->{'beta_barrel'} ne '') {
			print "<tr></td><td valign='top' align='left' colspan=29 bgcolor='$bg_colour'><span class='smaller'>Beta-barrel membrane segments are not counted as being membrane helical segments.<br>The inside/outside topology calculations do not include any beta-barrel membrane proteins or soluble proteins.</span></td></tr>\n";
		} elsif ($global->{'kingdom'}->{'soluble'} ne '') {
			print "<tr></td><td valign='top' align='left' colspan=29 bgcolor='$bg_colour'><span class='smaller'>The inside/outside topology calculations do not include any beta-barrel membrane proteins or soluble proteins.</span></td></tr>\n";
			print "<tr><td height=10 colspan=29 bgcolor='$bg_colour'></td></tr>\n";
		}
		print "</table>\n";

		print "</td>\n</tr>\n</table>\n";

		print "<br><br>\n";

	#===================================================================================================
	#	display the graph section
	#===================================================================================================

		print "<table cellspacing=0 cellpadding=10 border=1 bgcolor='#eeeeee'><tr><td valign='top' align='left'>\n";
		print "<table cellspacing=0 cellpadding=0 border=0>\n";
		print "<tr><td valign='top' align='left'><span class='h1'>GRAPH OF BENCHMARK RESULTS :</span><br><br>\n";

		print "<form action='benchmark.pl' method='post'>\n";
		print "<table border=0 cellpadding=0 cellspacing=0><tr><td>\n";

		print "<table border=0 cellpadding=0 cellspacing=0>\n";
		print "<tr><td align='left' valign='top'>\n";
		print "<input id='graph_coords' name='graph_coords' value='$graph_coords' type='hidden'>\n";
		print "<input id='graph_title' name='graph_title' value='$graph_title' type='hidden'>\n";
		print "<input id='graph_legend' name='graph_legend' value='$graph_legend' type='hidden'>\n";
		print "<input name='submit' value='Create Graph' type='submit'>\n";
		print "</td><td width=50></td>\n";
		print "<td align='left' valign='top'><span class='smaller'>&lt;== PLEASE NOTE : If you are viewing this from a Windows machine, your browser may not allow you to view this dynamically generated chart conveniently in the browser.<br>Instead, it may ask you to download the graph as a file and save it to disk.<br>When you do this, you may need to rename the ending file extension of the downloaded file to gif (eg. rename the file to benchmark.gif),<br>so that double-clicking this saved file on your disk will display this graph that has been dynamically generated for you.</span>\n";
		print "</td></tr></table>\n";

		print "</td></tr></table></form>\n";

		print "</td></tr></table></td></tr></table><br>\n"

	}

	#===================================================================================================
	#	dummy divs so that html body onload javascript will not fail
	#===================================================================================================

	my $output_line = "<div id='div_for_choose_standards'></div><div id='div_for_choose_standards_show'></div>\n";
	print $output_line;
	my $output_line = "<div id='div_for_more_dataset_options'></div><div id='div_for_more_dataset_options_show'></div>\n";
	print $output_line;
	my $output_line = "<div id='div_for_beta_barrel'></div><div id='div_for_beta_barrel_show'></div>\n";
	print $output_line;
	my $output_line = "<div id='div_for_soluble'></div><div id='div_for_soluble_show'></div>\n";
	print $output_line;
	my $output_line = "<div id='div_for_enter_seqs'></div><div id='div_for_enter_seqs_show'></div>\n";
	print $output_line;
	my $output_line = "<div id='div_for_upload_file'></div><div id='div_for_upload_file_show'></div>\n";
	print $output_line;
	my $output_line = "<div id='div_for_benchmark_criteria'></div><div id='div_for_benchmark_criteria_show'></div>\n";
	print $output_line;

	#===================================================================================================
	#	generic bottom of the displayed html page
	#===================================================================================================

	my $start_writing_html_from_file = 0;
	foreach my $input_line (@input_lines) {

		chomp $input_line;
		if ($input_line ne '') {

			if ($input_line eq "<a name='footer'>") {
				$start_writing_html_from_file = 1;
			}

			if ($input_line eq "<div id='div_for_key_to_benchmark_results_show'>") {
				$start_writing_html_from_file = 0;
			}
			elsif ($input_line eq "<a name='key_view'>") {
				$start_writing_html_from_file = 1;
			}
			elsif ($input_line eq "<div id='div_for_key_to_display_aligned_show'>") {
				$start_writing_html_from_file = 1;
			}

			if ($start_writing_html_from_file == 1) {
				my $output_line = $input_line . "\n";
				print $output_line;
			}
		}
	}

	if ($debug ne '') {
		print "<br>\n$debug<br>\n";
	}
}

#===================================================================================================
#	display_html_aligned_sequences
#===================================================================================================

sub display_html_aligned_sequences {

	print $q->header();

	#===================================================================================================
	#	generic top of the displayed html page
	#===================================================================================================

	my $html_file = 'benchmark.html';
	open HTMLFILE, $html_file or die $!;
	my @input_lines = <HTMLFILE>;

	my $stop_writing_html_from_file = 0;
	foreach my $input_line (@input_lines) {

		if ($stop_writing_html_from_file == 0) {

			chomp $input_line;
			if ($input_line ne '') {

				my $output_line = $input_line . "\n";
				print $output_line;

				if ($input_line eq "<a name='msg'>") {
					$stop_writing_html_from_file = 1;
				}
			}
		}
	}

	if ($global->{'info_msg'} ne '') {
		my $output_line = "<span class='red_bold'>" . $global->{'info_msg'} . "</span><br><br>\n\n";
		print $output_line;
	}

	#===================================================================================================
	#	display the relevent legend at the top of the page before showing the results
	#===================================================================================================

	my $start_writing_html_from_file = 0;
	foreach my $input_line (@input_lines) {

		chomp $input_line;
		if ($input_line ne '') {

			if ($input_line eq "<div id='div_for_key_to_display_aligned_show'>") {
				$start_writing_html_from_file = 1;
			}
			elsif ($input_line eq "<a name='end_key'>") {
				$start_writing_html_from_file = 0;
			}
			elsif ($input_line eq "<a name='benchmark_references'>") {
				$start_writing_html_from_file = 0;
			}

			if ($start_writing_html_from_file == 1) {
				my $output_line = $input_line . "\n";
				print $output_line;
			}
		}
	}
	my $output_line = "<br><br>\n\n";
	print $output_line;

	#===================================================================================================
	#	the displayed report view of aligned sequences and predictions
	#===================================================================================================

	#$debug .= Dumper( $output->{'display_method'} );
	#print $debug;

	my $output_line = "<span class='h1'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;DISPLAY OF ALIGNED SEQUENCES AND PREDICTIONS :</span><br><br><br>\n";
	print $output_line;

	foreach my $seqid ( sort keys %{$output->{'id'}->{'all'}} ) {

		my $pdbid = substr( $seqid, 0, 4 );
		my $aaseq = $output->{'aaseq'}->{$seqid};
		my $aaseq_length = length($aaseq);
		my $dssp = $output->{'dssp'}->{$seqid};
		my $name = $output->{'name'}->{$seqid};
		my $date = '';
		my $mpstruc_category = '';
		if ($global->{'run'} eq 'display_date') {
			$date = $output->{'release_dates'}->{$seqid};
			$mpstruc_category = $output->{'mpstruc_category'}->{$seqid};
		}
		my $blank_sequence = '_' x $aaseq_length;
		my $observed_helices = $blank_sequence;
		if (exists($output->{'observed_helices'}->{$seqid})) {
			$observed_helices = $output->{'observed_helices'}->{$seqid};
			if (($observed_helices eq '') or ($observed_helices eq '-')) {
				$observed_helices = $blank_sequence;
			}
		}
		my $opm_membrane_near = $blank_sequence;
		if (exists($output->{'opm_membrane_near'}->{$seqid})) {
			$opm_membrane_near = $output->{'opm_membrane_near'}->{$seqid};
			if (($opm_membrane_near eq '') or ($opm_membrane_near eq '-')) {
				$opm_membrane_near = $blank_sequence;
			}
		}
		my $opm_posn = $blank_sequence;
		if (exists($output->{'opm_posn'}->{$seqid})) {
			$opm_posn = $output->{'opm_posn'}->{$seqid};
			if (($opm_posn eq '') or ($opm_posn eq '-')) {
				$opm_posn = $blank_sequence;
			}
		}
		my $opm_seq = $output->{'opm_seq'}->{$seqid};
		my $opm_tm_segments = $output->{'opm_tm_segments'}->{$seqid};
		my $opm_tm_segments_topology = $blank_sequence;
		if (exists($output->{'opm_tm_segments_topology'}->{$seqid})) {
			$opm_tm_segments_topology = $output->{'opm_tm_segments_topology'}->{$seqid};
			if (($opm_tm_segments_topology eq '') or ($opm_tm_segments_topology eq '-')) {
				$opm_tm_segments_topology = $blank_sequence;
			}
		}
		my $opm_tmh_tmhh_segments = $blank_sequence;
		if (exists($output->{'opm_tmh_tmbb_segments'}->{$seqid})) {
			$opm_tmh_tmhh_segments = $output->{'opm_tmh_tmbb_segments'}->{$seqid};
			if (($opm_tmh_tmhh_segments eq '') or ($opm_tmh_tmhh_segments eq '-')) {
				$opm_tmh_tmhh_segments = $blank_sequence;
			}
		}
		my $opm_adj_tmh_tmhh_segments = $blank_sequence;
		if (exists($output->{'opm_adj_tmh_tmbb_segments'}->{$seqid})) {
			$opm_adj_tmh_tmhh_segments = $output->{'opm_adj_tmh_tmbb_segments'}->{$seqid};
			if (($opm_adj_tmh_tmhh_segments eq '') or ($opm_adj_tmh_tmhh_segments eq '-')) {
				$opm_adj_tmh_tmhh_segments = $blank_sequence;
			}
		}
		my $opm_tm_subunits = $blank_sequence;
		if (exists($output->{'opm_tm_subunits'}->{$seqid})) {
			$opm_tm_subunits = $output->{'opm_tm_subunits'}->{$seqid};
			if (($opm_tm_subunits eq '') or ($opm_tm_subunits eq '-')) {
				$opm_tm_subunits = $blank_sequence;
			}
		}
		my $opm_adj_tm_segments = $output->{'opm_adj_tm_segments'}->{$seqid};
		my $opm_adj_tm_segments_topology = $blank_sequence;
		if (exists($output->{'opm_adj_tm_segments_topology'}->{$seqid})) {
			$opm_adj_tm_segments_topology = $output->{'opm_adj_tm_segments_topology'}->{$seqid};
			if (($opm_adj_tm_segments_topology eq '') or ($opm_adj_tm_segments_topology eq '-')) {
				$opm_adj_tm_segments_topology = $blank_sequence;
			}
		}
		my $pdbtm_segments = $output->{'pdbtm_segments'}->{$seqid};
		my $pdbtm_segments_topology = $blank_sequence;
		if (exists($output->{'pdbtm_segments_topology'}->{$seqid})) {
			$pdbtm_segments_topology = $output->{'pdbtm_segments_topology'}->{$seqid};
			if (($pdbtm_segments_topology eq '') or ($pdbtm_segments_topology eq '-')) {
				$pdbtm_segments_topology = $blank_sequence;
			}
		}
		my $stride = $output->{'stride'}->{$seqid};
		my $stride_seq = $output->{'stride_seq'}->{$seqid};
		if (($stride eq '') or ($stride eq '-')) {
			$stride = $blank_sequence;
		}
		if (($stride_seq eq '') or ($stride_seq eq '-')) {
			$stride_seq = $blank_sequence;
		}

		my $output_line = $seqid . ' : ' . $name;
		if ($global->{'run'} eq 'display_date') {
			# $output_line .= " [$date] ($mpstruc_category)&nbsp;&nbsp;&nbsp;<a href='http://www.rcsb.org/pdb/explore/explore.do?structureId=$pdbid' target='_blank'>PDB link</a>&nbsp;&nbsp;&nbsp;<a href='http://opm.phar.umich.edu/protein.php?search=$pdbid' target='_blank'>OPM link</a>&nbsp;&nbsp;&nbsp;<a href='http://pdbtm.enzim.hu/?m=search&id=$pdbid' target='_blank'>PDBTM link</a>";
			my $this_pdb_link = 'http://www.rcsb.org/pdb/explore/explore.do?structureId=' . $pdbid;
			my $this_opm_link = 'http://opm.phar.umich.edu/protein.php?search=' . $pdbid;
			#my $this_pdbtm_link = 'http://pdbtm.enzim.hu/?_=/hitlist/Protein/%2560pdb_id%2560%2520LIKE%2520%2527%2525' . $pdbid . '%2525%2527&m=search&id=' . $pdbid;
			my $this_pdbtm_link = 'http://pdbtm_old.enzim.hu/?m=show_entry&id=' . $pdbid . '&nojava';
			$output_line .= " [$date] ($mpstruc_category)&nbsp;&nbsp;&nbsp;<a href='$this_pdb_link' target='_blank'>PDB link</a>&nbsp;&nbsp;&nbsp;<a href='$this_opm_link' target='_blank'>OPM link</a>&nbsp;&nbsp;&nbsp;<a href='$this_pdbtm_link' target='_blank'>PDBTM link</a>";
		}
		print $output_line . "<br>\n";
		my $sprintf_string = '%-' . $output->{'display_name_width'} . 's'; # "%-8s"
		my $display_width = $global->{'display_width'};
		if ($display_width == 0) {
			$display_width = $aaseq_length;
			if ($display_width == 0) {
				$display_width = 80;
			}
		}
		my $display_name_width = $output->{'display_name_width'};
		my $has_remainder = 0;
		if (($aaseq_length % $display_width) > 0) {
			$has_remainder = 1;
		}
		my $num_display_blocks_per_seq = int($aaseq_length / $display_width) + $has_remainder;

		foreach my $method (@{$output->{'display_method_list'}}) {
			if (exists($output->{'display_method'}->{$method})) {
				if ($output->{'display_method'}->{$method} == 1) {
					if (($output->{$method}->{$seqid} eq '') or ($output->{$method}->{$seqid} eq '-')) {
						$output->{$method}->{$seqid} = $blank_sequence;
					}
				}
			}
		}

		foreach my $method (@{$output->{'prediction_method'}}) {
			if (($output->{'prediction_sequence'}->{$method}->{$seqid} eq '') or ($output->{'prediction_sequence'}->{$method}->{$seqid} eq '-')) {
				$output->{'prediction_sequence'}->{$method}->{$seqid} = $blank_sequence;
			}
			if ( length($output->{'prediction_sequence'}->{$method}->{$seqid}) < $aaseq_length) {
				my $add_length = $aaseq_length - length($output->{'prediction_sequence'}->{$method}->{$seqid});
				my $add_string = '_' x $aaseq_length;
				$output->{'prediction_sequence'}->{$method}->{$seqid} .= $add_string;
			}
		}

		print "<span class='tt'>\n";

		for ( my $i = 0; $i < $num_display_blocks_per_seq; $i++ ) {

			my $display_from_residue = $i * $display_width + 1;
			my $display_to_residue = ($i + 1) * $display_width;
			if ($display_to_residue > $aaseq_length) {
				$display_to_residue = $aaseq_length;
			}
			my $substr_from_residue = $display_from_residue - 1;
			my $substr_length = $display_to_residue - $display_from_residue + 1;

			my $short_name;

			if ($global->{'dont_show_sequences'} == 0) {
				$short_name = substr( $data->{'method_short_name'}->{'pdb_seq'}, 0, $display_name_width );
				$short_name = sprintf( "%-${display_name_width}s", $short_name );
				$short_name =~ s/ /\&nbsp\;/g;
				$short_name = sprintf($sprintf_string, $short_name) . ':';
				$output_line = $short_name;
				$output_line .= substr( $aaseq, $substr_from_residue, $substr_length );
				print '<nobr>' . $output_line . "</nobr>&nbsp;<br>\n";
			}

			if (($global->{'show_intermediate'} == 1) and ($global->{'dont_show_sequences'} == 0)) {

				$short_name = substr( $data->{'method_short_name'}->{'opm_seq'}, 0, $display_name_width );
				$short_name = sprintf( "%-${display_name_width}s", $short_name );
				$short_name =~ s/ /\&nbsp\;/g;
				$short_name = sprintf($sprintf_string, $short_name) . ':';
				$output_line = $short_name;
				$output_line .= substr( $opm_seq, $substr_from_residue, $substr_length );
				print '<nobr>' . $output_line . "</nobr>&nbsp;<br>\n";

				$short_name = substr( $data->{'method_short_name'}->{'stride_seq'}, 0, $display_name_width );
				$short_name = sprintf( "%-${display_name_width}s", $short_name );
				$short_name =~ s/ /\&nbsp\;/g;
				$short_name = sprintf($sprintf_string, $short_name) . ':';
				$output_line = $short_name;
				$output_line .= substr( $stride_seq, $substr_from_residue, $substr_length );
				print '<nobr>' . $output_line . "</nobr>&nbsp;<br>\n";
			}

			$short_name = substr( 'BENCHMARK', 0, $display_name_width );
			$short_name = sprintf( "%-${display_name_width}s", $short_name );
			$short_name =~ s/ /\&nbsp\;/g;
			$short_name = sprintf($sprintf_string, $short_name) . ':';
			$output_line = $short_name;
			$output_line .= substr( $observed_helices, $substr_from_residue, $substr_length );
			print '<nobr>' . $output_line . "</nobr>&nbsp;<br>\n";

			foreach my $method (@{$output->{'display_method_list'}}) {
				if (exists($output->{'display_method'}->{$method})) {
					if ($output->{'display_method'}->{$method} == 1) {
						my $aaseq = $output->{$method}->{$seqid};
						$short_name = substr( $data->{'method_short_name'}->{$method}, 0, $display_name_width );
						$short_name = sprintf( "%-${display_name_width}s", $short_name );
						$short_name =~ s/ /\&nbsp\;/g;
						$short_name = sprintf($sprintf_string, $short_name) . ':';
						$output_line = $short_name;
						$output_line .= substr( $aaseq, $substr_from_residue, $substr_length );
						print '<nobr>' . $output_line . "</nobr>&nbsp;<br>\n";
					}
				}
			}

			if (($global->{'show_topology_files'} == 1) or ($global->{'show_intermediate'} == 1)) {

				$short_name = substr( $data->{'method_short_name'}->{'opm_adj_tm_segments_topology'}, 0, $display_name_width );
				$short_name = sprintf( "%-${display_name_width}s", $short_name );
				$short_name =~ s/ /\&nbsp\;/g;
				$short_name = sprintf($sprintf_string, $short_name) . ':';
				$output_line = $short_name;
				$output_line .= substr( $opm_adj_tm_segments_topology, $substr_from_residue, $substr_length );
				print '<nobr>' . $output_line . "</nobr>&nbsp;<br>\n";

				$short_name = substr( $data->{'method_short_name'}->{'opm_tm_segments_topology'}, 0, $display_name_width );
				$short_name = sprintf( "%-${display_name_width}s", $short_name );
				$short_name =~ s/ /\&nbsp\;/g;
				$short_name = sprintf($sprintf_string, $short_name) . ':';
				$output_line = $short_name;
				$output_line .= substr( $opm_tm_segments_topology, $substr_from_residue, $substr_length );
				print '<nobr>' . $output_line . "</nobr>&nbsp;<br>\n";

				if ($output->{'display_method'}->{'pdbtm_segments_topology'} == 0) { # if = 1 then it has already been displayed above
					$short_name = substr( $data->{'method_short_name'}->{'pdbtm_segments_topology'}, 0, $display_name_width );
					$short_name = sprintf( "%-${display_name_width}s", $short_name );
					$short_name =~ s/ /\&nbsp\;/g;
					$short_name = sprintf($sprintf_string, $short_name) . ':';
					$output_line = $short_name;
					$output_line .= substr( $pdbtm_segments_topology, $substr_from_residue, $substr_length );
					print '<nobr>' . $output_line . "</nobr>&nbsp;<br>\n";
				}
			}

			if ($global->{'show_intermediate'} == 1) {

				$short_name = substr( $data->{'method_short_name'}->{'opm_tm_subunits'}, 0, $display_name_width );
				$short_name = sprintf( "%-${display_name_width}s", $short_name );
				$short_name =~ s/ /\&nbsp\;/g;
				$short_name = sprintf($sprintf_string, $short_name) . ':';
				$output_line = $short_name;
				$output_line .= substr( $opm_tm_subunits, $substr_from_residue, $substr_length );
				print '<nobr>' . $output_line . "</nobr>&nbsp;<br>\n";

				$short_name = substr( $data->{'method_short_name'}->{'opm_tmh_tmbb_segments'}, 0, $display_name_width );
				$short_name = sprintf( "%-${display_name_width}s", $short_name );
				$short_name =~ s/ /\&nbsp\;/g;
				$short_name = sprintf($sprintf_string, $short_name) . ':';
				$output_line = $short_name;
				$output_line .= substr( $opm_tmh_tmhh_segments, $substr_from_residue, $substr_length );
				print '<nobr>' . $output_line . "</nobr>&nbsp;<br>\n";

				$short_name = substr( $data->{'method_short_name'}->{'opm_adj_tmh_tmbb_segments'}, 0, $display_name_width );
				$short_name = sprintf( "%-${display_name_width}s", $short_name );
				$short_name =~ s/ /\&nbsp\;/g;
				$short_name = sprintf($sprintf_string, $short_name) . ':';
				$output_line = $short_name;
				$output_line .= substr( $opm_adj_tmh_tmhh_segments, $substr_from_residue, $substr_length );
				print '<nobr>' . $output_line . "</nobr>&nbsp;<br>\n";

				$short_name = substr( $data->{'method_short_name'}->{'opm_posn'}, 0, $display_name_width );
				$short_name = sprintf( "%-${display_name_width}s", $short_name );
				$short_name =~ s/ /\&nbsp\;/g;
				$short_name = sprintf($sprintf_string, $short_name) . ':';
				$output_line = $short_name;
				$output_line .= substr( $opm_posn, $substr_from_residue, $substr_length );
				print '<nobr>' . $output_line . "</nobr>&nbsp;<br>\n";

				$short_name = substr( $data->{'method_short_name'}->{'opm_membrane_near'}, 0, $display_name_width );
				$short_name = sprintf( "%-${display_name_width}s", $short_name );
				$short_name =~ s/ /\&nbsp\;/g;
				$short_name = sprintf($sprintf_string, $short_name) . ':';
				$output_line = $short_name;
				$output_line .= substr( $opm_membrane_near, $substr_from_residue, $substr_length );
				print '<nobr>' . $output_line . "</nobr>&nbsp;<br>\n";
			}

			foreach my $method (@{$output->{'prediction_method'}}) {

				my $prediction_sequence = $output->{'prediction_sequence'}->{$method}->{$seqid};
				$short_name = substr( $data->{'method_short_name'}->{$method}, 0, $display_name_width );
				$short_name = sprintf( "%-${display_name_width}s", $short_name );
				$short_name =~ s/ /\&nbsp\;/g;
				$short_name = sprintf($sprintf_string, $short_name) . ':';
				$output_line = $short_name;
				$output_line .= substr( $prediction_sequence, $substr_from_residue, $substr_length );
				print '<nobr>' . $output_line . "</nobr>&nbsp;<br>\n";
			}
			print "<br>\n";
		}
		print "</span>\n";
		print "<br><br>\n";
	}

	#===================================================================================================
	#	dummy divs so that html body onload javascript will not fail
	#===================================================================================================

	my $output_line = "<div id='div_for_choose_standards'></div><div id='div_for_choose_standards_show'></div>\n";
	print $output_line;
	my $output_line = "<div id='div_for_more_dataset_options'></div><div id='div_for_more_dataset_options_show'></div>\n";
	print $output_line;
	my $output_line = "<div id='div_for_beta_barrel'></div><div id='div_for_beta_barrel_show'></div>\n";
	print $output_line;
	my $output_line = "<div id='div_for_soluble'></div><div id='div_for_soluble_show'></div>\n";
	print $output_line;
	my $output_line = "<div id='div_for_enter_seqs'></div><div id='div_for_enter_seqs_show'></div>\n";
	print $output_line;
	my $output_line = "<div id='div_for_upload_file'></div><div id='div_for_upload_file_show'></div>\n";
	print $output_line;
	my $output_line = "<div id='div_for_benchmark_criteria'></div><div id='div_for_benchmark_criteria_show'></div>\n";
	print $output_line;

	#===================================================================================================
	#	generic bottom of the displayed html page - without the legend that was already displayed above
	#===================================================================================================

	my $start_writing_html_from_file = 0;
	foreach my $input_line (@input_lines) {

		chomp $input_line;
		if ($input_line ne '') {

			if ($input_line eq "<a name='footer'>") {
				$start_writing_html_from_file = 1;
			}

			elsif ($input_line eq "<div id='div_for_key_to_display_aligned_show'>") {
				$start_writing_html_from_file = 0;
			}
			elsif ($input_line eq "<a name='end_key'>") {
				$start_writing_html_from_file = 1;
			}
			elsif ($input_line eq "<a name='benchmark_references'>") {
				$start_writing_html_from_file = 1;
			}

			if ($start_writing_html_from_file == 1) {
				my $output_line = $input_line . "\n";
				print $output_line;
			}
		}
	}

	if ($debug ne '') {
		print "<br>\n$debug<br>\n";
	}
}

#===================================================================================================
#	build_list_of_output_ids
#===================================================================================================

sub build_list_of_output_ids {

	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel', 'soluble' );
	my @protein_categories2 = ( 'alpha_polytopic_bitopic', 'beta_barrel' );

	my %output_id;
	$output->{'id'} = \%output_id;
	my %output_id_all;
	my %output_id_alpha_polytopic_bitopic;
	my %output_id_beta_barrel;
	my %output_id_soluble;
	$output->{'id'}->{'all'} = \%output_id_all;
	$output->{'id'}->{'alpha_polytopic_bitopic'} = \%output_id_alpha_polytopic_bitopic;
	$output->{'id'}->{'beta_barrel'} = \%output_id_beta_barrel;
	$output->{'id'}->{'soluble'} = \%output_id_soluble;

	my %output_flag1; # contains all possible ids (for this protein category), and flag equals 1 if this id has survived all criteria
	my %output_flag2; # contains all possible ids (for this protein category), and flag equals 0 ready for recording next set of criteria
	$output->{'flag1'} = \%output_flag1;
	$output->{'flag2'} = \%output_flag2;
	my %output_flag1_alpha_polytopic_bitopic;
	my %output_flag1_beta_barrel;
	my %output_flag1_soluble;
	my %output_flag2_alpha_polytopic_bitopic;
	my %output_flag2_beta_barrel;
	my %output_flag2_soluble;
	$output->{'flag1'}->{'alpha_polytopic_bitopic'} = \%output_flag1_alpha_polytopic_bitopic;
	$output->{'flag1'}->{'beta_barrel'} = \%output_flag1_beta_barrel;
	$output->{'flag1'}->{'soluble'} = \%output_flag1_soluble;
	$output->{'flag2'}->{'alpha_polytopic_bitopic'} = \%output_flag2_alpha_polytopic_bitopic;
	$output->{'flag2'}->{'beta_barrel'} = \%output_flag2_beta_barrel;
	$output->{'flag2'}->{'soluble'} = \%output_flag2_soluble;

	#===================================================================================================
	# this is the list of all protein-chain ids
	#===================================================================================================

	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {
			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'id'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					$output->{'flag1'}->{$this_protein_category}->{$this_id} = 1;
					$output->{'flag2'}->{$this_protein_category}->{$this_id} = 0;
				}
			}
		}
	}

	#===================================================================================================
	# restrict output list by similarity criteria
	#===================================================================================================

	for my $this_protein_category ( @protein_categories2 ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {
			if ($global->{'nonhomologous_list'}->{$this_protein_category} ne '') {
				my $file_having_list_of_ids = $global->{'nonhomologous_list'}->{$this_protein_category};
				open LISTFILE, $file_having_list_of_ids or die $!;
				my @input_lines = <LISTFILE>;
				foreach my $input_line (@input_lines) {
					chomp $input_line;
					if ($input_line ne '') {
						my @bits = split(/:::::/, $input_line);
						my $this_id = trim($bits[0]);
						$output->{'flag2'}->{$this_protein_category}->{$this_id} = 1;
					}
				}
				foreach my $this_id ( keys %{$output->{'flag1'}->{$this_protein_category}} ) {
					if (($output->{'flag1'}->{$this_protein_category}->{$this_id} == 1) && ($output->{'flag2'}->{$this_protein_category}->{$this_id} == 1)) {
						$output->{'flag1'}->{$this_protein_category}->{$this_id} = 1;
					} else {
						$output->{'flag1'}->{$this_protein_category}->{$this_id} = 0;
					}
					$output->{'flag2'}->{$this_protein_category}->{$this_id} = 0 ;
				}
			}
		}
	}

	#===================================================================================================
	# restrict output list by kingdom criteria
	#===================================================================================================

	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {
			my $search_kingdom = '';
			if ($global->{'kingdom'}->{$this_protein_category} eq 'eukaryote') {
				$search_kingdom = 'Superkingdom: Eukaryota';
			} elsif ($global->{'kingdom'}->{$this_protein_category} eq 'bacteria') {
				$search_kingdom = 'Superkingdom: Bacteria';
			} elsif ($global->{'kingdom'}->{$this_protein_category} eq 'archaea') {
				$search_kingdom = 'Superkingdom: Archaea';
			} elsif ($global->{'kingdom'}->{$this_protein_category} eq 'viruses') {
				$search_kingdom = 'Superkingdom: Viruses';
			}
			if ($search_kingdom ne '') {
				my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'pdbe_taxonomy'};
				open LISTFILE, $file_having_list_of_ids or die $!;
				my @input_lines = <LISTFILE>;
				foreach my $input_line (@input_lines) {
					chomp $input_line;
					if ($input_line ne '') {
						my @bits = split(/:::::/, $input_line);
						my $this_id = trim($bits[0]);
						my $this_taxonomy_list = trim($bits[1]);
						if ( $this_taxonomy_list =~ /$search_kingdom/ ) {
							$output->{'flag2'}->{$this_protein_category}->{$this_id} = 1;
						}
					}
				}
				foreach my $this_id ( keys %{$output->{'flag1'}->{$this_protein_category}} ) {
					if (($output->{'flag1'}->{$this_protein_category}->{$this_id} == 1) && ($output->{'flag2'}->{$this_protein_category}->{$this_id} == 1)) {
						$output->{'flag1'}->{$this_protein_category}->{$this_id} = 1;
					} else {
						$output->{'flag1'}->{$this_protein_category}->{$this_id} = 0;
					}
					$output->{'flag2'}->{$this_protein_category}->{$this_id} = 0 ;
				}
			}
		}
	}

	#===================================================================================================
	# restrict output list by polytopic/monotopic criteria
	#===================================================================================================

	if ($global->{'kingdom'}->{'alpha_polytopic_bitopic'} ne '') {
		my $this_protein_category = 'alpha_polytopic_bitopic';
		if ($global->{'tmh_numTMH'} ne '') {
			my $search_polytopic_monotopic = '';
			if ($global->{'tmh_numTMH'} eq 'polytopic') {
				$search_polytopic_monotopic = 'OPM-POLYTOPIC';
			} elsif ($global->{'tmh_numTMH'} eq 'bitopic') {
				$search_polytopic_monotopic = 'OPM-BITOPIC';
			} elsif ($global->{'tmh_numTMH'} eq 'polytopic_bitopic') {
				$search_polytopic_monotopic = '(OPM-POLYTOPIC)|(OPM-BITOPIC)';
			}
			if ($search_polytopic_monotopic ne '') {
				my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'opm_category'};
				open LISTFILE, $file_having_list_of_ids or die $!;
				my @input_lines = <LISTFILE>;
				foreach my $input_line (@input_lines) {
					chomp $input_line;
					if ($input_line ne '') {
						my @bits = split(/:::::/, $input_line);
						my $this_id = trim($bits[0]);
						my $this_polytopic_monotopic = trim($bits[1]);
						if ( $this_polytopic_monotopic =~ /$search_polytopic_monotopic/ ) {
							$output->{'flag2'}->{$this_protein_category}->{$this_id} = 1;
						}
					}
				}
				foreach my $this_id ( keys %{$output->{'flag1'}->{$this_protein_category}} ) {
					if (($output->{'flag1'}->{$this_protein_category}->{$this_id} == 1) && ($output->{'flag2'}->{$this_protein_category}->{$this_id} == 1)) {
						$output->{'flag1'}->{$this_protein_category}->{$this_id} = 1;
					} else {
						$output->{'flag1'}->{$this_protein_category}->{$this_id} = 0;
					}
					$output->{'flag2'}->{$this_protein_category}->{$this_id} = 0 ;
				}
			}
		}
	}

	#===================================================================================================
	# restrict output list by helix type criteria
	#===================================================================================================

	if ($global->{'kingdom'}->{'alpha_polytopic_bitopic'} ne '') {
		my $this_protein_category = 'alpha_polytopic_bitopic';
		if ($global->{'helix_type'} ne '') {
			my $search_helix_type = '';
			if ($global->{'helix_type'} eq 'nohalfmembrane') {
				$search_helix_type = 'nohalfmembrane';
			} elsif ($global->{'helix_type'} eq 'halfmembrane') {
				$search_helix_type = 'halfmembrane';
			}
			if ($search_helix_type ne '') {
				foreach my $this_id ( keys %{$output->{'flag1'}->{$this_protein_category}} ) {
					$output->{'flag2'}->{$this_protein_category}->{$this_id} = 0;
				}
				my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'reentrant_helices'};
				open LISTFILE, $file_having_list_of_ids or die $!;
				my @input_lines = <LISTFILE>;
				foreach my $input_line (@input_lines) {
					chomp $input_line;
					if ($input_line ne '') {
						my @bits = split(/:::::/, $input_line);
						my $this_id = trim($bits[0]);
						if ($this_id ne '-') {
							$output->{'flag2'}->{$this_protein_category}->{$this_id} = 1;
						}
					}
				}
				foreach my $this_id ( keys %{$output->{'flag1'}->{$this_protein_category}} ) {
					if ($search_helix_type eq 'halfmembrane') {
						if (($output->{'flag1'}->{$this_protein_category}->{$this_id} == 1) && ($output->{'flag2'}->{$this_protein_category}->{$this_id} == 1)) {
							$output->{'flag1'}->{$this_protein_category}->{$this_id} = 1;
						} else {
							$output->{'flag1'}->{$this_protein_category}->{$this_id} = 0;
						}
					} elsif ($search_helix_type eq 'nohalfmembrane') {
						if (($output->{'flag1'}->{$this_protein_category}->{$this_id} == 1) && ($output->{'flag2'}->{$this_protein_category}->{$this_id} == 0)) {
							$output->{'flag1'}->{$this_protein_category}->{$this_id} = 1;
						} else {
							$output->{'flag1'}->{$this_protein_category}->{$this_id} = 0;
						}
					}
					$output->{'flag2'}->{$this_protein_category}->{$this_id} = 0 ;
				}
			}
			if ($global->{'helix_type'} eq 'none') {
				foreach my $this_id ( keys %{$output->{'flag1'}->{$this_protein_category}} ) {
					$output->{'flag1'}->{$this_protein_category}->{$this_id} = 0;
				}
			}
		}
	}

	#===================================================================================================
	# restrict output list by family criteria
	#===================================================================================================

	if ($global->{'kingdom'}->{'alpha_polytopic_bitopic'} ne '') {
		my $this_protein_category = 'alpha_polytopic_bitopic';
		if ($global->{'tmh_family'} ne '') {
			my $search_family = $global->{'tmh_family'};
			$search_family =~ s/mpstruc__//;
			if ($search_family eq 'all') {
				$global->{'tmh_family'} = '';
			} else {
				$search_family =~ s/\(/\\\(/;
				$search_family =~ s/\)/\\\)/;
				$search_family =~ s/\//\\\//;
				$search_family =~ s/\+/\\\+/;
				$search_family =~ s/\-/\\\-/;
			}
			if ($global->{'tmh_family'} ne '') {
				my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'mpstruc_category'};
				open LISTFILE, $file_having_list_of_ids or die $!;
				my @input_lines = <LISTFILE>;
				foreach my $input_line (@input_lines) {
					chomp $input_line;
					if ($input_line ne '') {
						my @bits = split(/:::::/, $input_line);
						my $this_id = trim($bits[0]);
						my $this_polytopic_monotopic = trim($bits[1]);
						$this_polytopic_monotopic =~ s/\(/\\\(/;
						$this_polytopic_monotopic =~ s/\)/\\\)/;
						$this_polytopic_monotopic =~ s/\//\\\//;
						$this_polytopic_monotopic =~ s/\+/\\\+/;
						$this_polytopic_monotopic =~ s/\-/\\\-/;
						#if ( $this_polytopic_monotopic =~ /$search_family/ ) {
						if ( $this_polytopic_monotopic eq $search_family ) {
							$output->{'flag2'}->{$this_protein_category}->{$this_id} = 1;
						}
					}
				}
				foreach my $this_id ( keys %{$output->{'flag1'}->{$this_protein_category}} ) {
					if (($output->{'flag1'}->{$this_protein_category}->{$this_id} == 1) && ($output->{'flag2'}->{$this_protein_category}->{$this_id} == 1)) {
						$output->{'flag1'}->{$this_protein_category}->{$this_id} = 1;
					} else {
						$output->{'flag1'}->{$this_protein_category}->{$this_id} = 0;
					}
					$output->{'flag2'}->{$this_protein_category}->{$this_id} = 0 ;
				}
			}
		}
	}

	#===================================================================================================
	# restrict output list by PDB homologue dates
	#===================================================================================================

	if ($global->{'kingdom'}->{'alpha_polytopic_bitopic'} ne '') {
		if ($global->{'restrict_year'} ne '') {
			my $this_protein_category = 'alpha_polytopic_bitopic';
			my $search_year = $global->{'restrict_year'} . '-01-01';
			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'earliest_release_date'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					my $this_homologue_date = $bits[1];
					if ($this_homologue_date lt $global->{'restrict_year'}) {
						$output->{'flag2'}->{$this_protein_category}->{$this_id} = 0;
					} else {
						$output->{'flag2'}->{$this_protein_category}->{$this_id} = 1;
					}
				}
			}
			foreach my $this_id ( keys %{$output->{'flag1'}->{$this_protein_category}} ) {
				if (($output->{'flag1'}->{$this_protein_category}->{$this_id} == 1) && ($output->{'flag2'}->{$this_protein_category}->{$this_id} == 1)) {
					$output->{'flag1'}->{$this_protein_category}->{$this_id} = 1;
				} else {
					$output->{'flag1'}->{$this_protein_category}->{$this_id} = 0;
				}
				$output->{'flag2'}->{$this_protein_category}->{$this_id} = 0 ;
			}
		}
	}

	#===================================================================================================
	# restrict output list by resolution or experimental method
	#===================================================================================================

	if (($global->{'restrict_experimental_method'} ne '') or ($global->{'restrict_resolution'} ne '')) {
		for my $this_protein_category ( @protein_categories ) {
			if ($global->{'kingdom'}->{$this_protein_category} ne '') {
				my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'pdb_method_resolution'};
				open LISTFILE, $file_having_list_of_ids or die $!;
				my @input_lines = <LISTFILE>;
				foreach my $input_line (@input_lines) {
					chomp $input_line;
					if ($input_line ne '') {
						# 1mhs_A:::::METHOD,,,,,ELECTRON CRYSTALLOGRAPHY;;;;;RESOLUTION,,,,,8.0;;;;;:::::
						# 1n7l_A:::::METHOD,,,,,SOLUTION NMR;;;;;:::::
						my @bits = split(/:::::/, $input_line);
						my $this_id = trim($bits[0]);
						my $this_experimental_method = '';
						my $this_resolution = '';
						if ( exists ($output->{'flag1'}->{$this_protein_category}->{$this_id}) ) {
							if ($output->{'flag1'}->{$this_protein_category}->{$this_id} == 1) {
								my @bits2 = split(/\;\;\;\;\;/, $bits[1]);
								foreach my $this_method_resolution (@bits2) {
									my @bits2 = split(/\,\,\,\,\,/, $this_method_resolution);
									my $this_name = $bits2[0];
									my $this_value = $bits2[1];
									if ($this_name eq 'METHOD') {
										$this_experimental_method = $this_value;
									} elsif ($this_name = 'RESOLUTION') {
										$this_resolution = $this_value;
									}
								}
							}
						}
						my $flag2_1 = 0;
						my $flag2_2 = 0;
						if ($global->{'restrict_experimental_method'} eq '') {
							$flag2_1 = 1;
						} else {
							my $search_experimental_method = ';;;;;' . $this_experimental_method . ';;;;;';
							if ( index( $global->{'restrict_experimental_method'}, $search_experimental_method ) > -1 ) {
								$flag2_1 = 1;
							}
						}
						if ($global->{'restrict_resolution'} eq '') {
							$flag2_2 = 1;
						} else {
							if ($this_resolution eq '') {
								$flag2_2 = 1;
							} else {
								if ($this_resolution <= $global->{'restrict_resolution'}) {
									$flag2_2 = 1;
								}
							}
						}
						if (($flag2_1 == 1) && ($flag2_2 == 1)) {
							$output->{'flag2'}->{$this_protein_category}->{$this_id} = 1;
						}
					}
				}
				foreach my $this_id ( keys %{$output->{'flag1'}->{$this_protein_category}} ) {
					if (($output->{'flag1'}->{$this_protein_category}->{$this_id} == 1) && ($output->{'flag2'}->{$this_protein_category}->{$this_id} == 1)) {
						$output->{'flag1'}->{$this_protein_category}->{$this_id} = 1;
					} else {
						$output->{'flag1'}->{$this_protein_category}->{$this_id} = 0;
					}
					$output->{'flag2'}->{$this_protein_category}->{$this_id} = 0 ;
				}
			}
		}
	}

	#===================================================================================================
	# restrict the number of soluble proteins in the list if a number was given
	#===================================================================================================

	my $this_protein_category = 'soluble';
	if ($global->{'kingdom'}->{$this_protein_category} ne '') {
		if (($global->{'num_soluble'} ne 'all') and ($global->{'num_soluble'} ne '')) {
			my $count_number_in_list = 0;
			my $max_number_in_list = $global->{'num_soluble'};
			foreach my $this_id ( keys %{$output->{'flag1'}->{$this_protein_category}} ) {
				if ($output->{'flag1'}->{$this_protein_category}->{$this_id} == 1) {
					if ($count_number_in_list < $max_number_in_list) {
						$count_number_in_list++;
					} else {
						$output->{'flag1'}->{$this_protein_category}->{$this_id} = 0;
					}
				}
			}
		}
	}

	#===================================================================================================
	# finalise the list
	#===================================================================================================

	for my $this_protein_category ( @protein_categories ) {
		foreach my $this_id ( keys %{$output->{'flag1'}->{$this_protein_category}} ) {
			if ($output->{'flag1'}->{$this_protein_category}->{$this_id} == 1) {
				$output->{'id'}->{$this_protein_category}->{$this_id} = '';
				$output->{'id'}->{'all'}->{$this_id} = '';
			}
		}
	}
	# $debug .= Dumper( $output->{'id'}->{'all'} );
}

#===================================================================================================
#	build_text_results
#===================================================================================================

sub build_text_results {

	#===================================================================================================
	# user has requested to see FASTA sequences for the requested set of protein-chain ids
	#===================================================================================================

	if (($global->{'run'} eq 'fasta') or ($global->{'run'} eq 'fastaG') or ($global->{'run'} eq 'dload') or ($global->{'run'} eq 'dloadG')) {

		get_aaseq_for_output();
	}

	#===================================================================================================
	# user has requested to see name (and release-date) for the requested set of protein-chain ids
	#===================================================================================================

	if ($global->{'run'} eq 'names') {

		get_names_for_output();

		get_release_dates_for_output();

		get_methods_resolutions_for_output();

		get_mpstruc_category_for_output();
	}
}

#===================================================================================================
#	build_benchmark_results
#===================================================================================================

sub build_benchmark_results {

	my %output_id = %{$output->{'id'}->{'all'}};
	foreach my $this_id ( keys %output_id ) {
		$output->{'num_reference_seqs'} += 1;
	}

	# read in the sequences to be benchmarked

	get_aaseq_for_output();

	# now that we have the sequences, validate that user entered predictions are not too long

	if ($global->{'use_the_user_predictions'} == 1) {
		validate_user_prediction_sequence_lengths();
	}

	# read in the benchmark data

	if ($global->{'benchmark_standard'} eq 'opm') {
		get_opm_tm_segments_for_output( 'benchmark' );
		$output->{'observed_helices'} = $output->{'opm_tm_segments'};
		$output->{'observed_topology'} = $output->{'opm_tm_segments_topology'};
	}
	if (($global->{'benchmark_standard'} eq 'pdbtm_helices_loops') or ($global->{'benchmark_standard'} eq 'pdbtm_helices')) {
		get_pdbtm_segments_for_output( 'benchmark' );
		$output->{'observed_helices'} = $output->{'pdbtm_segments'};
		$output->{'observed_topology'} = $output->{'pdbtm_segments_topology'};
	}
	if ($global->{'benchmark_standard'} eq 'opm_adj') {
		get_opm_adj_tm_segments_for_output( 'benchmark' );
		$output->{'observed_helices'} = $output->{'opm_adj_tm_segments'};
		$output->{'observed_topology'} = $output->{'opm_adj_tm_segments_topology'};
	}
	$output->{'observed_segments'} = convert_sequences_to_membrane_segments( $output->{'observed_helices'}, $global->{'benchmark_standard'} );

	foreach my $this_id ( keys %{$output->{'observed_segments'}} ) {
		my $this_segments = $output->{'observed_segments'}->{$this_id};
		if ($this_segments ne '-') {
			my @bits = split(/,/, $this_segments);
			$output->{'num_reference_membrane_helices'} += @bits;
		}
	}

	# for one of the helix-end scores, calculate a calculation-intensive score only once

	$global->{'helix_ends_score_denominator'} = 1;
	$global->{'PI'} = 3.1415926539;
	my $variance = $global->{'max_dist_approx_ends'};
	my $variance2 = 2 * $variance;
	my $std_dev = sqrt $variance;
	my $std_dev_root_2_pi = $std_dev * (sqrt (2 * $global->{'PI'}));
	$global->{'helix_ends_score_denominator'} = 1 / $std_dev_root_2_pi;
	$global->{'std_dev'} = $std_dev;
	$global->{'variance'} = $variance;
	$global->{'variance2'} = $variance2;
	$global->{'std_dev_root_2_pi'} = $std_dev_root_2_pi;

	# read in the prediction data for each prediction method
	# then calculate the benchmark statistics for it

	foreach my $method (@{$output->{'prediction_method'}}) {

		my %output_prediction_sequences;

		if ($method eq 'you') {

			foreach my $this_id ( keys %output_id ) {
				$output_prediction_sequences{$this_id} = $data->{'user_data'}->{$this_id};
			}

		} else {

			my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel', 'soluble' );
			for my $this_protein_category ( @protein_categories ) {

				my $input_file = $data->{'file_name'}->{$this_protein_category}->{$method};
				open INFILE, $input_file or die $!;
				my @input_lines = <INFILE>;
				foreach my $input_line (@input_lines) {
					chomp $input_line;
					if ($input_line ne '') {
						my @bits = split(/:::::/, $input_line);
						my $this_id = $bits[0];
						my $this_prediction_sequence = $bits[1];

						if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {

							#if ($this_prediction_sequence eq '-') {
							#	$this_prediction_sequence = '_' x length( $output->{'aaseq'}->{$this_id} );
							#}

							$output_prediction_sequences{$this_id} = $this_prediction_sequence;
						}
					}
				}
			}
		}

		my $output_prediction_segments = convert_sequences_to_membrane_segments( \%output_prediction_sequences, $method );

		my $stats = calculate_benchmark_stats( $output->{'aaseq'}, $output->{'observed_segments'}, $output->{'observed_topology'}, $output_prediction_segments, \%output_prediction_sequences, $method );

		push( @{$output->{'qok'}}, $stats->{'qok'} );
		push( @{$output->{'qhtm_obs'}}, $stats->{'qhtm_obs'} );
		push( @{$output->{'qhtm_prd'}}, $stats->{'qhtm_prd'} );
		push( @{$output->{'avhe_diff'}}, $stats->{'avhe_diff'} );
		push( @{$output->{'qhe_obs'}}, $stats->{'qhe_obs'} );
		push( @{$output->{'qhe_prd'}}, $stats->{'qhe_prd'} );
		push( @{$output->{'qnhe_obs'}}, $stats->{'qnhe_obs'} );
		push( @{$output->{'qnhe_prd'}}, $stats->{'qnhe_prd'} );
		push( @{$output->{'q2'}}, $stats->{'q2'} );
		push( @{$output->{'htm_mcc'}}, $stats->{'htm_mcc'} );
		push( @{$output->{'q2t_obs'}}, $stats->{'q2t_obs'} );
		push( @{$output->{'q2t_prd'}}, $stats->{'q2t_prd'} );
		push( @{$output->{'q2n_obs'}}, $stats->{'q2n_obs'} );
		push( @{$output->{'q2n_prd'}}, $stats->{'q2n_prd'} );
		push( @{$output->{'qok3'}}, $stats->{'qok3'} );
		push( @{$output->{'nterm'}}, $stats->{'nterm'} );
		push( @{$output->{'q2_ioseg'}}, $stats->{'q2_ioseg'} );
		push( @{$output->{'qiom_obs'}}, $stats->{'qiom_obs'} );
		push( @{$output->{'qiom_prd'}}, $stats->{'qiom_prd'} );
		push( @{$output->{'qio_obs'}}, $stats->{'qio_obs'} );
		push( @{$output->{'qio_prd'}}, $stats->{'qio_prd'} );
		push( @{$output->{'q3'}}, $stats->{'q3'} );
		push( @{$output->{'q2_iores'}}, $stats->{'q2_iores'} );
		push( @{$output->{'io_mcc'}}, $stats->{'io_mcc'} );
	}

	# sort the output results by one of the display columns

	my @output_methods = @{$output->{'prediction_method'}};
	my $sort_by = trim( $global->{'sort_by'} );
	if ($sort_by ne '') {

		my @output_qok = @{$output->{'qok'}};
		my @output_qhtm_obs = @{$output->{'qhtm_obs'}};
		my @output_qhtm_prd = @{$output->{'qhtm_prd'}};
		my @output_avhe_diff = @{$output->{'avhe_diff'}};
		my @output_qhe_obs = @{$output->{'qhe_obs'}};
		my @output_qhe_prd = @{$output->{'qhe_prd'}};
		my @output_qnhe_obs = @{$output->{'qnhe_obs'}};
		my @output_qnhe_prd = @{$output->{'qnhe_prd'}};
		my @output_q2 = @{$output->{'q2'}};
		my @output_htm_mcc = @{$output->{'htm_mcc'}};
		my @output_q2t_obs = @{$output->{'q2t_obs'}};
		my @output_q2t_prd = @{$output->{'q2t_prd'}};
		my @output_q2n_obs = @{$output->{'q2n_obs'}};
		my @output_q2n_prd = @{$output->{'q2n_prd'}};
		my @output_qok3 = @{$output->{'qok3'}};
		my @output_nterm = @{$output->{'nterm'}};
		my @output_q2_ioseg = @{$output->{'q2_ioseg'}};
		my @output_qiom_obs = @{$output->{'qiom_obs'}};
		my @output_qiom_prd = @{$output->{'qiom_prd'}};
		my @output_qio_obs = @{$output->{'qio_obs'}};
		my @output_qio_prd = @{$output->{'qio_prd'}};
		my @output_q3 = @{$output->{'q3'}};
		my @output_q2_iores = @{$output->{'q2_iores'}};
		my @output_io_mcc = @{$output->{'io_mcc'}};

		my @sort_lines;
		my $qok_i = 0;
		my $qhtm_obs_i = 3;
		my $qhtm_prd_i = 6;
		my $qhe_obs_i = 9;
		my $qhe_prd_i = 12;
		my $qnhe_obs_i = 15;
		my $qnhe_prd_i = 18;
		my $q2_i = 21;
		my $q2t_obs_i = 24;
		my $q2t_prd_i = 27;
		my $q2n_obs_i = 30;
		my $q2n_prd_i = 33;
		my $qok3_i = 36;
		my $nterm_i = 39;
		my $q2_ioseg_i = 42;
		my $qiom_obs_i = 45;
		my $qiom_prd_i = 48;
		my $qio_obs_i = 51;
		my $qio_prd_i = 54;
		my $q3_i = 57;
		my $q2_iores_i = 60;
		my $htm_mcc_i = 63;
		my $io_mcc_i = 76;
		my $avhe_diff_i = 89;
		my $method_i = 102;

		my $high_number_for_inverting_sort = 999;

		for ( my $i = 0; $i < @output_methods; $i++ ) {
			my $method = sprintf("%-13s", $output_methods[$i] );
			my $qok = sprintf("%03d", $output_qok[$i] );
			my $qhtm_obs = sprintf("%03d", $output_qhtm_obs[$i] );
			my $qhtm_prd = sprintf("%03d", $output_qhtm_prd[$i] );
			my $qhe_obs = sprintf("%03d", $output_qhe_obs[$i] );
			my $qhe_prd = sprintf("%03d", $output_qhe_prd[$i] );
			my $qnhe_obs = sprintf("%03d", $output_qnhe_obs[$i] );
			my $qnhe_prd = sprintf("%03d", $output_qnhe_prd[$i] );
			my $q2 = sprintf("%03d", $output_q2[$i] );
			my $q2t_obs = sprintf("%03d", $output_q2t_obs[$i] );
			my $q2t_prd = sprintf("%03d", $output_q2t_prd[$i] );
			my $q2n_obs = sprintf("%03d", $output_q2n_obs[$i] );
			my $q2n_prd = sprintf("%03d", $output_q2n_prd[$i] );
			my $qok3 = sprintf("%03d", $output_qok3[$i] );
			my $nterm = sprintf("%03d", $output_nterm[$i] );
			my $q2_ioseg = sprintf("%03d", $output_q2_ioseg[$i] );
			my $qiom_obs = sprintf("%03d", $output_qiom_obs[$i] );
			my $qiom_prd = sprintf("%03d", $output_qiom_prd[$i] );
			my $qio_obs = sprintf("%03d", $output_qio_obs[$i] );
			my $qio_prd = sprintf("%03d", $output_qio_prd[$i] );
			my $q3 = sprintf("%03d", $output_q3[$i] );
			my $q2_iores = sprintf("%03d", $output_q2_iores[$i] );

			my $extra_padding = '';
			my $htm_mcc = sprintf("%2.10f", $output_htm_mcc[$i] );
			for ( my $j = 0; $j < (13 - length($htm_mcc)); $j++ ) {
				$extra_padding .= '0';
			}
			$htm_mcc = $extra_padding . $htm_mcc;

			$extra_padding = '';
			my $io_mcc = sprintf("%2.10f", $output_io_mcc[$i] );
			for ( my $j = 0; $j < (13 - length($io_mcc)); $j++ ) {
				$extra_padding .= '0';
			}
			$io_mcc = $extra_padding . $io_mcc;

			$extra_padding = '';
			my $avhe_diff = $high_number_for_inverting_sort - $output_avhe_diff[$i];
			$avhe_diff = sprintf("%3.9f", $avhe_diff );
			for ( my $j = 0; $j < (13 - length($avhe_diff)); $j++ ) {
				$extra_padding .= '0';
			}
			$avhe_diff = $extra_padding . $avhe_diff;

			# default : $sort_by eq 'qok'
			my $sort_line = $qok . $qhtm_obs . $qhtm_prd . $qhe_obs . $qhe_prd . $qnhe_obs . $qnhe_prd . $q2 . $q2t_obs . $q2t_prd . $q2n_obs . $q2n_prd . $qok3 . $nterm . $q2_ioseg . $qiom_obs . $qiom_prd . $qio_obs . $qio_prd . $q3 . $q2_iores . $htm_mcc . $io_mcc . $avhe_diff . $method;

			if ($sort_by eq 'qhtm_obs') {
				$sort_line = $qhtm_obs . '__________' . $sort_line;
			} elsif ($sort_by eq 'qhtm_prd') {
				$sort_line = $qhtm_prd . '__________' . $sort_line;
			} elsif ($sort_by eq 'qhe_obs') {
				$sort_line = $qhe_obs . '__________' . $sort_line;
			} elsif ($sort_by eq 'qhe_prd') {
				$sort_line = $qhe_prd . '__________' . $sort_line;
			} elsif ($sort_by eq 'qnhe_obs') {
				$sort_line = $qnhe_obs . '__________' . $sort_line;
			} elsif ($sort_by eq 'qnhe_prd') {
				$sort_line = $qnhe_prd . '__________' . $sort_line;
			} elsif ($sort_by eq 'q2') {
				$sort_line = $q2 . '__________' . $sort_line;
			} elsif ($sort_by eq 'htm_mcc') {
				$sort_line = $htm_mcc . $sort_line;
			} elsif ($sort_by eq 'q2t_obs') {
				$sort_line = $q2t_obs . '__________' . $sort_line;
			} elsif ($sort_by eq 'q2t_prd') {
				$sort_line = $q2t_prd . '__________' . $sort_line;
			} elsif ($sort_by eq 'q2n_obs') {
				$sort_line = $q2n_obs . '__________' . $sort_line;
			} elsif ($sort_by eq 'q2n_prd') {
				$sort_line = $q2n_prd . '__________' . $sort_line;
			} elsif ($sort_by eq 'qok3') {
				$sort_line = $qok3 . '__________' . $sort_line;
			} elsif ($sort_by eq 'nterm') {
				$sort_line = $nterm . '__________' . $sort_line;
			} elsif ($sort_by eq 'q2_ioseg') {
				$sort_line = $q2_ioseg . '__________' . $sort_line;
			} elsif ($sort_by eq 'qiom_obs') {
				$sort_line = $qiom_obs . '__________' . $sort_line;
			} elsif ($sort_by eq 'qiom_prd') {
				$sort_line = $qiom_prd . '__________' . $sort_line;
			} elsif ($sort_by eq 'qio_obs') {
				$sort_line = $qio_obs . '__________' . $sort_line;
			} elsif ($sort_by eq 'qio_prd') {
				$sort_line = $qio_prd . '__________' . $sort_line;
			} elsif ($sort_by eq 'q3') {
				$sort_line = $q3 . '__________' . $sort_line;
			} elsif ($sort_by eq 'q2_iores') {
				$sort_line = $q2_iores . '__________' . $sort_line;
			} elsif ($sort_by eq 'io_mcc') {
				$sort_line = $io_mcc . $sort_line;
			} elsif ($sort_by eq 'avhe_diff') {
				$sort_line = $avhe_diff . $sort_line;
			} elsif ($sort_by eq 'alphabet') {
				my $this_method = $output_methods[$i];
				my $padded_method = substr( $data->{'method_long_name'}->{$this_method}, 0, 53 );
				$padded_method = uc( $padded_method );
				$padded_method = sprintf("%-53s", $padded_method );
				$sort_line = $padded_method . $sort_line;
			} else {
				$sort_line = '___' . '__________' . $sort_line;
			}
			push( @sort_lines, $sort_line );
		}

		#$debug .= "sort_by=$sort_by <br><br>\n\n";#####
		#foreach my $debug_sort (@sort_lines) {#####
		#	$debug .= "BEFORE : $debug_sort<BR>\n";#####
		#}#####

		my @sorted_lines;
		if ($sort_by eq 'alphabet') {
			@sorted_lines = sort {$a cmp $b} @sort_lines;
		} else {
			# reverse/descending sort
			@sorted_lines = sort {$b cmp $a} @sort_lines;
		}
		#$debug .= Dumper( @sorted_lines ) . "<br>";

		#foreach my $debug_sort (@sorted_lines) {#####
		#	$debug .= "AFTER : $debug_sort<BR>\n";#####
		#}#####

		{
			my @output_method;
			my @output_qok;
			my @output_qhtm_obs;
			my @output_qhtm_prd;
			my @output_avhe_diff;
			my @output_qhe_obs;
			my @output_qhe_prd;
			my @output_qnhe_obs;
			my @output_qnhe_prd;
			my @output_q2;
			my @output_htm_mcc;
			my @output_q2t_obs;
			my @output_q2t_prd;
			my @output_q2n_obs;
			my @output_q2n_prd;
			my @output_qok3;
			my @output_nterm;
			my @output_q2_ioseg;
			my @output_qiom_obs;
			my @output_qiom_prd;
			my @output_qio_obs;
			my @output_qio_prd;
			my @output_q3;
			my @output_q2_iores;
			my @output_io_mcc;
			$output->{'method'} = \@output_method;
			$output->{'qok'} = \@output_qok;
			$output->{'qhtm_obs'} = \@output_qhtm_obs;
			$output->{'qhtm_prd'} = \@output_qhtm_prd;
			$output->{'avhe_diff'} = \@output_avhe_diff;
			$output->{'qhe_obs'} = \@output_qhe_obs;
			$output->{'qhe_prd'} = \@output_qhe_prd;
			$output->{'qnhe_obs'} = \@output_qnhe_obs;
			$output->{'qnhe_prd'} = \@output_qnhe_prd;
			$output->{'q2'} = \@output_q2;
			$output->{'htm_mcc'} = \@output_htm_mcc;
			$output->{'q2t_obs'} = \@output_q2t_obs;
			$output->{'q2t_prd'} = \@output_q2t_prd;
			$output->{'q2n_obs'} = \@output_q2n_obs;
			$output->{'q2n_prd'} = \@output_q2n_prd;
			$output->{'qok3'} = \@output_qok3;
			$output->{'nterm'} = \@output_nterm;
			$output->{'q2_ioseg'} = \@output_q2_ioseg;
			$output->{'qiom_obs'} = \@output_qiom_obs;
			$output->{'qiom_prd'} = \@output_qiom_prd;
			$output->{'qio_obs'} = \@output_qio_obs;
			$output->{'qio_prd'} = \@output_qio_prd;
			$output->{'q3'} = \@output_q3;
			$output->{'q2_iores'} = \@output_q2_iores;
			$output->{'io_mcc'} = \@output_io_mcc;

			foreach my $sorted_line (@sorted_lines) {
				#if ($sort_by eq 'alphabet') {
				#	$sorted_line = substr( $sorted_line, 100 );
				#} else {
				#	$sorted_line = substr( $sorted_line, 13 );
				#}
				if ($sort_by eq 'alphabet') {
					$sorted_line = substr( $sorted_line, 53 );
				} else {
					$sorted_line = substr( $sorted_line, 13 );
				}

				my $method = substr( $sorted_line, $method_i );
				$method = trim($method);
				my $method_name = $output->{'method_name'}->{$method};

				my $qok = substr( $sorted_line, $qok_i, 3 );
				my $qhtm_obs = substr( $sorted_line, $qhtm_obs_i, 3 );
				my $qhtm_prd = substr( $sorted_line, $qhtm_prd_i, 3 );
				my $qhe_prd = substr( $sorted_line, $qhe_prd_i, 3 );
				my $qnhe_obs = substr( $sorted_line, $qnhe_obs_i, 3 );
				my $qnhe_prd = substr( $sorted_line, $qnhe_prd_i, 3 );
				my $q2 = substr( $sorted_line, $q2_i, 3 );
				my $htm_mcc = substr( $sorted_line, $htm_mcc_i, 13 );
				my $q2t_obs = substr( $sorted_line, $q2t_obs_i, 3 );
				my $q2t_prd = substr( $sorted_line, $q2t_prd_i, 3 );
				my $q2n_obs = substr( $sorted_line, $q2n_obs_i, 3 );
				my $q2n_prd = substr( $sorted_line, $q2n_prd_i, 3 );
				my $qhe_obs = substr( $sorted_line, $qhe_obs_i, 3 );
				my $qok3 = substr( $sorted_line, $qok3_i, 3 );
				my $nterm = substr( $sorted_line, $nterm_i, 3 );
				my $q2_ioseg = substr( $sorted_line, $q2_ioseg_i, 3 );
				my $qiom_obs = substr( $sorted_line, $qiom_obs_i, 3 );
				my $qiom_prd = substr( $sorted_line, $qiom_prd_i, 3 );
				my $qio_obs = substr( $sorted_line, $qio_obs_i, 3 );
				my $qio_prd = substr( $sorted_line, $qio_prd_i, 3 );
				my $q3 = substr( $sorted_line, $q3_i, 3 );
				my $q2_iores = substr( $sorted_line, $q2_iores_i, 3 );
				my $io_mcc = substr( $sorted_line, $io_mcc_i, 13 );
				my $avhe_diff = $high_number_for_inverting_sort - substr( $sorted_line, $avhe_diff_i, 13 );

				$qok = sprintf("%d", $qok );
				$qhtm_obs = sprintf("%d", $qhtm_obs );
				$qhtm_prd = sprintf("%d", $qhtm_prd );
				$avhe_diff = sprintf("%.1f", $avhe_diff );
				$qhe_obs = sprintf("%d", $qhe_obs );
				$qhe_prd = sprintf("%d", $qhe_prd );
				$qnhe_obs = sprintf("%d", $qnhe_obs );
				$qnhe_prd = sprintf("%d", $qnhe_prd );
				$q2 = sprintf("%d", $q2 );
				$htm_mcc = sprintf("%.2f", $htm_mcc );
				$q2t_obs = sprintf("%d", $q2t_obs );
				$q2t_prd = sprintf("%d", $q2t_prd );
				$q2n_obs = sprintf("%d", $q2n_obs );
				$q2n_prd = sprintf("%d", $q2n_prd );

				my $this_method_has_topology = 0;
				foreach my $topology_method ( @{$data->{'list_of_topology_prediction_methods'}} ) {
					if ($method eq $topology_method) {
						$this_method_has_topology = 1;
					}
				}
				if ($this_method_has_topology == 1) {
					$qok3 = sprintf("%d", $qok3 );
					$nterm = sprintf("%d", $nterm );
					$q2_ioseg = sprintf("%d", $q2_ioseg );
					$qiom_obs = sprintf("%d", $qiom_obs );
					$qiom_prd = sprintf("%d", $qiom_prd );
					$qio_obs = sprintf("%d", $qio_obs );
					$qio_prd = sprintf("%d", $qio_prd );
					$q3 = sprintf("%d", $q3 );
					$q2_iores = sprintf("%d", $q2_iores );
					$io_mcc = sprintf("%.2f", $io_mcc );
				} else {
					$qok3 = '';
					$nterm = '';
					$q2_ioseg = '';
					$qiom_obs = '';
					$qiom_prd = '';
					$qio_obs = '';
					$qio_prd = '';
					$q3 = '';
					$q2_iores = '';
					$io_mcc = '';
				}

				push( @output_method, $method );
				push( @output_qok, $qok );
				push( @output_qhtm_obs, $qhtm_obs );
				push( @output_qhtm_prd, $qhtm_prd );
				push( @output_avhe_diff, $avhe_diff );
				push( @output_qhe_obs, $qhe_obs );
				push( @output_qhe_prd, $qhe_prd );
				push( @output_qnhe_obs, $qnhe_obs );
				push( @output_qnhe_prd, $qnhe_prd );
				push( @output_q2, $q2 );
				push( @output_htm_mcc, $htm_mcc );
				push( @output_q2t_obs, $q2t_obs );
				push( @output_q2t_prd, $q2t_prd );
				push( @output_q2n_obs, $q2n_obs );
				push( @output_q2n_prd, $q2n_prd );
				push( @output_qok3, $qok3 );
				push( @output_nterm, $nterm );
				push( @output_q2_ioseg, $q2_ioseg );
				push( @output_qiom_obs, $qiom_obs );
				push( @output_qiom_prd, $qiom_prd );
				push( @output_qio_obs, $qio_obs );
				push( @output_qio_prd, $qio_prd );
				push( @output_q3, $q3 );
				push( @output_q2_iores, $q2_iores );
				push( @output_io_mcc, $io_mcc );
			}

			$output->{'prediction_method'} = \@output_method;
			#$debug .= Dumper( @output_method ) . "<br>";
		}
	}
}

#===================================================================================================
#	convert_sequences_to_membrane_segments
#===================================================================================================

sub convert_sequences_to_membrane_segments {

	my $aaseq_hash = shift;
	my %aaseq = %$aaseq_hash;
	my $this_method = shift;
	my $return_segments;
	my %return_segments_hash;
	$return_segments = \%return_segments_hash;

	# pdbtm_helices_loops, pdbtm_helices	:	2a0l_A:::::ALPHA__side_1,,,,,52-156,225-237,;;;;;ALPHA__side_2,,,,,24-29,174-178,206-206,;;;;;ALPHA__alpha_helix,,,,,30-51,157-173,207-224,;;;;;ALPHA__membrane_loop,,,,,179-205,;;;;;:::::

	foreach my $this_id ( keys %{$aaseq_hash} ) {

		my $this_sequence = $aaseq_hash->{$this_id};

		my $this_return_segments = '';

		$this_sequence =~ s/H/M/g; # membrane helix
		$this_sequence =~ s/X/M/g; # membrane helix in PDBTM, MEMSAT3
		$this_sequence =~ s/T/M/g; # membrane segment in PHDThtm_PBIL
		$this_sequence =~ s/R/M/g; # reentrant loop
		$this_sequence =~ s/L/M/g; # reentrant loop in PDBTM
		#$this_sequence =~ s/S/\_/g; # signal peptide

		# opm, opm_adj	:	2l34_A:::::oooooooMMMMMMMMMMMMMMMMMMMMMMimii:::::
		my @chars = split(//, $this_sequence);
		my $prev_char = '_';
		my $i = 1;
		my $this_from = 0;
		my $this_to = 0;
		foreach my $this_char (@chars) {
			if (($this_char eq 'M') and ($prev_char ne 'M')) {
				$this_from = $i;
			} elsif (($this_char ne 'M') and ($prev_char eq 'M')) {
				$this_to = $i - 1;
				$this_return_segments .= $this_from . '-' . $this_to . ',';
			}
			$prev_char = $this_char;
			$i++;
		}

		$return_segments->{$this_id} = $this_return_segments;
	}

	return $return_segments;
}

#===================================================================================================
#	build_view_results
#===================================================================================================

sub build_view_results {

	my %output_id = %{$output->{'id'}->{'all'}};
	my @output_display_method_list;
	$output->{'display_method_list'} = \@output_display_method_list;
	my %output_display_method = %{$output->{'display_method'}};

	$output->{'display_name_width'} = 9;

	get_aaseq_for_output();

	get_names_for_output();

	if ($global->{'run'} eq 'display_date') {

		get_release_dates_for_output();

		get_mpstruc_category_for_output();
	}

	#$debug .= Dumper( $global->{'benchmark_standard'} );

	if ($global->{'benchmark_standard'} eq 'opm_adj') {
		push( @output_display_method_list, 'opm_adj_tm_segments' );
	} elsif ($global->{'benchmark_standard'} eq 'opm') {
		push( @output_display_method_list, 'opm_tm_segments' );
	} elsif (($global->{'benchmark_standard'} eq 'pdbtm_helices_loops') or ($global->{'benchmark_standard'} eq 'pdbtm_helices')) {
		push( @output_display_method_list, 'pdbtm_segments' );
	}
	#-----
	if (($output->{'display_method'}->{'opm_adj_tm_segments'} == 1) or ($output->{'display_method'}->{'opm_adj_tm_segments_topology'} == 1) or ($global->{'benchmark_standard'} eq 'opm_adj') or ($global->{'show_topology_files'} == 1) or ($global->{'show_intermediate'} == 1)) {
		get_opm_adj_tm_segments_for_output( 'view' );
		if ($global->{'benchmark_standard'} ne 'opm_adj') {
			push( @output_display_method_list, 'opm_adj_tm_segments' );
		}
		if ($global->{'benchmark_standard'} eq 'opm_adj') {
			$output->{'observed_helices'} = $output->{'opm_adj_tm_segments'};
		}
		if (($output->{'display_method'}->{'opm_adj_tm_segments_topology'} == 1) or ($global->{'show_topology_files'} == 1) or ($global->{'show_intermediate'} == 1)) {
			push( @output_display_method_list, 'opm_adj_tm_segments_topology' );
		}
	}
	#-----
	if (($output->{'display_method'}->{'opm_tm_segments'} == 1) or ($output->{'display_method'}->{'opm_tm_segments_topology'} == 1) or ($global->{'benchmark_standard'} eq 'opm') or ($global->{'show_topology_files'} == 1) or ($global->{'show_intermediate'} == 1)) {
		get_opm_tm_segments_for_output( 'view' );
		if ($global->{'benchmark_standard'} ne 'opm') {
			push( @output_display_method_list, 'opm_tm_segments' );
		}
		if ($global->{'benchmark_standard'} eq 'opm') {
			$output->{'observed_helices'} = $output->{'opm_tm_segments'};
		}
		if (($output->{'display_method'}->{'opm_tm_segments_topology'} == 1) or ($global->{'show_topology_files'} == 1) or ($global->{'show_intermediate'} == 1)) {
			push( @output_display_method_list, 'opm_tm_segments_topology' );
		}
	}
	#-----
	if (($output->{'display_method'}->{'pdbtm_segments'} == 1) or ($output->{'display_method'}->{'pdbtm_segments_topology'} == 1) or ($global->{'benchmark_standard'} eq 'pdbtm_helices_loops') or ($global->{'benchmark_standard'} eq 'pdbtm_helices') or ($global->{'show_topology_files'} == 1) or ($global->{'show_intermediate'} == 1)) {
		get_pdbtm_segments_for_output( 'view' );
		if (($global->{'benchmark_standard'} ne 'pdbtm_helices_loops') and ($global->{'benchmark_standard'} ne 'pdbtm_helices')) {
			push( @output_display_method_list, 'pdbtm_segments' );
		}
		if (($global->{'benchmark_standard'} eq 'pdbtm_helices_loops') or ($global->{'benchmark_standard'} eq 'pdbtm_helices')) {
			$output->{'observed_helices'} = $output->{'pdbtm_segments'};
		}
		if (($output->{'display_method'}->{'pdbtm_segments_topology'} == 1) or ($global->{'show_topology_files'} == 1) or ($global->{'show_intermediate'} == 1)) {
			push( @output_display_method_list, 'pdbtm_segments_topology' );
		}
	}
	#-----
	if ($output->{'display_method'}->{'dssp'} == 1) {
		get_PDB_DSSP_for_output();
		push( @output_display_method_list, 'dssp' );
	}
	if ($output->{'display_method'}->{'stride'} == 1) {
		get_PDB_STRIDE_for_output();
		push( @output_display_method_list, 'stride' );
	}
	#-----
	if ($global->{'show_intermediate'} == 1) {
		get_PDB_file_seq_for_output();
		get_PDB_OPM_STRIDE_file_seq_for_output();
		get_OPM_tm_subunits_for_output();
		get_OPM_tmh_tmbb_for_output();
		get_OPM_adj_tmh_tmbb_for_output();
		get_OPM_topology_for_output();
		get_OPM_membrane_near_for_output();
	}
	#-----
	if (($global->{'show_topology_files'} == 1) or ($global->{'show_intermediate'} == 1)) {
		get_OPM_tm_segments_topology_for_output();
		get_OPM_adj_tm_segments_topology_for_output();
	}

	foreach my $method (@{$output->{'prediction_method'}}) {

		my %output_prediction_sequence_method_id;

		if ($method eq 'you') {

			foreach my $this_id ( keys %output_id ) {
				$output_prediction_sequence_method_id{$this_id} = $data->{'user_data'}->{$this_id};
			}

		} else {

			my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel', 'soluble' );
			for my $this_protein_category ( @protein_categories ) {

				my $get_method_for_this_protein_category = 0;
				if ($global->{'manual_list_seqs'} ne '') {
					$get_method_for_this_protein_category = 1;
				} else {
					if ($global->{'kingdom'}->{$this_protein_category} ne '') {
						$get_method_for_this_protein_category = 1;
					}
				}

				if ($get_method_for_this_protein_category == 1) {
					my $input_file = $data->{'file_name'}->{$this_protein_category}->{$method};
					open INFILE, $input_file or die $!;
					my @input_lines = <INFILE>;
					foreach my $input_line (@input_lines) {
						chomp $input_line;
						if ($input_line ne '') {
							my @bits = split(/:::::/, $input_line);
							my $this_id = $bits[0];
							my $this_prediction_sequence = $bits[1];

							if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {

								if ($this_prediction_sequence eq '-') {
									$this_prediction_sequence = '_' x length( $output->{'aaseq'}->{$this_id} );
								}

								$output_prediction_sequence_method_id{$this_id} = $this_prediction_sequence;
							}
						}
					}
				}
			}
		}

		$output->{'prediction_sequence'}->{$method} = \%output_prediction_sequence_method_id;
	}
}

#===================================================================================================
#	calculate_benchmark_stats
#===================================================================================================

sub calculate_benchmark_stats {

	my $reference_aaseq_hash = shift;
	my $reference_segments_hash = shift;
	my $reference_topology_hash = shift;
	my $method_segments_hash = shift;
	my $method_topology_hash = shift;
	my $method = shift;
	my %reference_aaseq = %$reference_aaseq_hash;
	my %reference_segments = %$reference_segments_hash;
	my %reference_topology = %$reference_topology_hash;
	my %method_segments = %$method_segments_hash;
	my %method_topology = %$method_topology_hash;

	my $stats;
	my %stats_hash;
	$stats = \%stats_hash;

	# initialise the displayed statistics totalled for all sequences in the prediction method

	$stats->{'qok'} = 0;
	$stats->{'qhtm_obs'} = 0;
	$stats->{'qhtm_prd'} = 0;
	$stats->{'avhe_diff'} = 0;
	$stats->{'qhe_obs'} = 0;
	$stats->{'qhe_prd'} = 0;
	$stats->{'qnhe_obs'} = 0;
	$stats->{'qnhe_prd'} = 0;
	$stats->{'q2'} = 0;
	$stats->{'htm_mcc'} = 0;
	$stats->{'q2t_obs'} = 0;
	$stats->{'q2t_prd'} = 0;
	$stats->{'q2n_obs'} = 0;
	$stats->{'q2n_prd'} = 0;
	$stats->{'qok3'} = 0;
	$stats->{'nterm'} = 0;
	$stats->{'q2_ioseg'} = 0;
	$stats->{'qiom_obs'} = 0;
	$stats->{'qiom_prd'} = 0;
	$stats->{'qio_obs'} = 0;
	$stats->{'qio_prd'} = 0;
	$stats->{'q3'} = 0;
	$stats->{'q2_iores'} = 0;
	$stats->{'io_mcc'} = 0;

	# initialise the overall counts that will be used to calculate the displayed statistics

	my $stats_num_seqs = 0;
	my $stats_num_obs_tms_in_this_method = 0;
	my $stats_num_prd_tms_in_this_method = 0;
	my $stats_num_correct_obs_tms_in_this_method = 0;
	my $stats_num_correct_prd_tms_in_this_method = 0;
	my $stats_num_correct_obs_and_prd_tms_in_this_method = 0;
	my $stats_num_obs_tms_residues_in_this_method = 0;
	my $stats_num_prd_tms_residues_in_this_method = 0;
	my $stats_num_obs_nontms_residues_in_this_method = 0;
	my $stats_num_prd_nontms_residues_in_this_method = 0;
	my $stats_num_correct_obs_tms_residues_in_this_method = 0;
	my $stats_num_correct_prd_tms_residues_in_this_method = 0;
	my $stats_num_correct_obs_nontms_residues_in_this_method = 0;
	my $stats_num_correct_prd_nontms_residues_in_this_method = 0;
	my $stats_ratio_correct_residues_in_this_method = 0;
	my $stats_sum_dist_obs_helix_end_minus_prd_helix_end = 0;
	my $stats_num_obs_approx_ends_in_this_method = 0;
	my $stats_num_prd_approx_ends_in_this_method = 0;
	my $stats_num_correct_obs_approx_ends_in_this_method = 0;
	my $stats_num_correct_prd_approx_ends_in_this_method = 0;
	my $stats_num_obs_scaled_ends_in_this_method = 0;
	my $stats_num_prd_scaled_ends_in_this_method = 0;
	my $stats_sum_of_scaled_scores_of_helix_ends = 0;
	my $stats_num_obs_iosegment_in_this_method = 0;
	my $stats_num_prd_iosegment_in_this_method = 0;
	my $stats_num_correct_obs_iosegment_in_this_method = 0;
	my $stats_num_correct_prd_iosegment_in_this_method = 0;
	my $stats_num_correct_obs_and_prd_iosegment_in_this_method = 0;
	my $stats_num_obs_iomsegment_in_this_method = 0;
	my $stats_num_prd_iomsegment_in_this_method = 0;
	my $stats_num_correct_obs_iomsegment_in_this_method = 0;
	my $stats_num_correct_prd_iomsegment_in_this_method = 0;
	my $stats_ratio_correct_iotopology_residues_in_this_method = 0;
	my $stats_ratio_correct_iomtopology_residues_in_this_method = 0;
	my $stats_ratio_correct_iotopology_segments_in_this_method = 0;
	my $stats_num_correct_nterminal_topology = 0;
	my $stats_qok3_count_for_this_method = 0;
	my $stats_num_i_topology_residues_in_this_method = 0;
	my $stats_num_o_topology_residues_in_this_method = 0;
	my $stats_num_correct_prd_i_topology_residues_in_this_method = 0;
	my $stats_num_correct_prd_o_topology_residues_in_this_method = 0;

	# does this method have topology (inside/outside) predictions, or only membrane segment (membrane helices) predictions

	my $this_method_has_topology = 0;
	foreach my $topology_method ( @{$data->{'list_of_topology_prediction_methods'}} ) {
		if ($method eq $topology_method) {
			$this_method_has_topology = 1;
		}
	}

	# for each sequence in the reference dataset, calculate the benchmark statistics

	foreach my $this_seqid (keys %method_segments) {

		my $this_aaseq = $reference_aaseq{$this_seqid};
		my $this_reference_segments = $reference_segments{$this_seqid};
		my $this_method_segments = $method_segments{$this_seqid};
		my $this_reference_topology = $reference_topology{$this_seqid};
		my $this_method_topology = $method_topology{$this_seqid};
		my $stats_num_residues = length($this_aaseq);
		if ($this_aaseq eq '-') {
			$stats_num_residues = 0;
		}
		my $stats_num_known_topology_residues = 0;
		my $stats_num_known_io_topology_residues = 0;

		$stats_num_seqs = $stats_num_seqs + 1;

		my $stats_qok_for_this_seq = 0;
		my $stats_num_obs_tms_in_this_seq = 0;
		my $stats_num_prd_tms_in_this_seq = 0;
		my $stats_num_correct_obs_tms_in_this_seq = 0;
		my $stats_num_correct_prd_tms_in_this_seq = 0;
		my $stats_dist_obs_helix_end_minus_prd_helix_end = 0;
		my $stats_qok3_for_this_seq = 0;
		my $stats_num_obs_iosegment_in_this_seq = 0;
		my $stats_num_prd_iosegment_in_this_seq = 0;
		my $stats_num_prd_iosegment_in_this_seq_with_interface_as_nonmembrane = 0;
		my $stats_num_correct_obs_iosegment_in_this_seq = 0;
		my $stats_num_correct_prd_iosegment_in_this_seq = 0;

		# list the observed helices for this sequence, observed by the benchmark reference
		my @observed_helices;
		my @observed_helices_start;
		my @observed_helices_end;
		my @observed_helices_already_counted;
		my @observed_residues;
		for ( my $i = 1; $i <= $stats_num_residues; $i++ ) {
			$observed_residues[$i] = 0;
		}
		if ($this_reference_segments ne '-') {
			@observed_helices = split(/\,/, $this_reference_segments);
			$stats_num_obs_tms_in_this_seq = @observed_helices;
			$stats_num_obs_tms_in_this_method += $stats_num_obs_tms_in_this_seq;
			$stats_num_obs_approx_ends_in_this_method += ($stats_num_obs_tms_in_this_seq * 2);
			$stats_num_obs_scaled_ends_in_this_method += ($stats_num_obs_tms_in_this_seq * 2);
			for ( my $i = 0; $i < @observed_helices; $i++ ) {
				my $observed_helix = $observed_helices[$i];
				$observed_helix =~ /(-?\d+)-(-?\d+)/;
				my $start = $1;
				my $end = $2;
				$observed_helices_start[$i] = $start;
				$observed_helices_end[$i] = $end;
				$observed_helices_already_counted[$i] = 0;
				for ( my $j = $start; $j <= $end; $j++ ) {
					$observed_residues[$j] = 1;
				}
			}
		}

		# list the observed topologies for this sequence, observed by the benchmark reference
		my @observed_topology_start;
		my @observed_topology_end;
		my @observed_topology_location;
		my @observed_topology_already_counted;
		my @observed_topology_residues;
		my $this_protein_has_topology = 0;
		if ($this_method_has_topology == 1) {
			for ( my $i = 1; $i <= $stats_num_residues; $i++ ) {
				$observed_topology_residues[$i] = '_';
			}
			if (($this_reference_topology ne '-') and ($this_reference_topology ne '')) {
				my @shifted_observed_topology_residues = split(//, $this_reference_topology);
				for ( my $i = 1; $i <= $stats_num_residues; $i++ ) {
					$observed_topology_residues[$i] = $shifted_observed_topology_residues[ $i - 1 ];
					if (($observed_topology_residues[$i] eq 'M') or ($observed_topology_residues[$i] eq 'H')) {
						$this_protein_has_topology = 1;
					}
				}
				if ($this_protein_has_topology == 1) {
					my $prev_char = '';
					my $prev_start = 0;
					my $prev_end = 0;
					for ( my $i = 1; $i <= $stats_num_residues; $i++ ) {
						my $this_char = $observed_topology_residues[$i];
						if ($this_char ne '_') {
							$stats_num_known_topology_residues++;
							if ($this_char eq 'i') {
								$stats_num_i_topology_residues_in_this_method++;
								$stats_num_known_io_topology_residues++;
							} elsif ($this_char eq 'o') {
								$stats_num_o_topology_residues_in_this_method++;
								$stats_num_known_io_topology_residues++;
							}
						}
						if ($this_char ne $prev_char) {
							# process the previous segment
							if ( ($prev_start != 0) and ($prev_char ne '_') and (($prev_char eq 'i') or ($prev_char eq 'o')) ) {
								$prev_end = $i - 1;
								push( @observed_topology_start, $prev_start );
								push( @observed_topology_end, $prev_end );
								push( @observed_topology_location, $prev_char );
								push( @observed_topology_already_counted, 0 );
							}
							# start processing of this new segment
							if (($this_char eq 'i') or ($this_char eq 'o')) {
								$prev_start = $i;
							}
						}
						$prev_char = $this_char;
					}
					# process the last segment
					if ( ($prev_start != 0) and ($prev_char ne '_') and (($prev_char eq 'i') or ($prev_char eq 'o')) ) {
						$prev_end = $stats_num_residues;
						push( @observed_topology_start, $prev_start );
						push( @observed_topology_end, $prev_end );
						push( @observed_topology_location, $prev_char );
						push( @observed_topology_already_counted, 0 );
					}
					$stats_num_obs_iosegment_in_this_seq = @observed_topology_start;
					$stats_num_obs_iosegment_in_this_method += $stats_num_obs_iosegment_in_this_seq;
				}
			}
		}

		# list the predicted helices for this sequence, predicted by the method to be benchmarked
		my @predicted_helices;
		my @predicted_helices_start;
		my @predicted_helices_end;
		my @predicted_helices_already_counted;
		my @predicted_residues;
		for ( my $i = 1; $i <= $stats_num_residues; $i++ ) {
			$predicted_residues[$i] = 0;
		}
		if (($this_method_segments ne '-') and ($this_method_segments ne '')) {
			@predicted_helices = split(/\,/, $this_method_segments);
			$stats_num_prd_tms_in_this_seq = @predicted_helices;
			$stats_num_prd_tms_in_this_method += $stats_num_prd_tms_in_this_seq;
			$stats_num_prd_approx_ends_in_this_method += ($stats_num_prd_tms_in_this_seq * 2);
			$stats_num_prd_scaled_ends_in_this_method += ($stats_num_prd_tms_in_this_seq * 2);
			for ( my $i = 0; $i < @predicted_helices; $i++ ) {
				my $predicted_helix = $predicted_helices[$i];
				$predicted_helix =~ /(-?\d+)-(-?\d+)/;
				my $start = $1;
				my $end = $2;
				$predicted_helices_start[$i] = $start;
				$predicted_helices_end[$i] = $end;
				$predicted_helices_already_counted[$i] = 0;
				for ( my $j = $start; $j <= $end; $j++ ) {
					$predicted_residues[$j] = 1;
				}
			}
		}

		# list the predicted topologies for this sequence, predicted by the method to be benchmarked
		my @predicted_topology_start;
		my @predicted_topology_end;
		my @predicted_topology_location;
		my @predicted_topology_residues;
		my @predicted_topology_start_with_interface_as_nonmembrane;
		my @predicted_topology_end_with_interface_as_nonmembrane;
		my @predicted_topology_location_with_interface_as_nonmembrane;
		my @predicted_topology_residues_with_interface_as_nonmembrane;
		if (($this_method_has_topology == 1) and ($this_protein_has_topology == 1)) {

			for ( my $i = 1; $i <= $stats_num_residues; $i++ ) {
				$predicted_topology_residues[$i] = '_';
				$predicted_topology_residues_with_interface_as_nonmembrane[$i] = '_';
			}
			my $this_method_topology_with_interface_as_nonmembrane = $this_method_topology; 	# used only for octopus and memsat3.
														# they predict interface helices that will be counted as membrane or non-membrane,
														# whichever makes the prediction look correct, ie whichever is the same as observed.

			if (($this_method_topology ne '-') and ($this_method_topology ne '')) {
				if (($method eq 'memsatsvm') or ($method eq 'philius') or ($method eq 'phobius') or ($method eq 'polyphobius') or ($method eq 'svmtop') or ($method eq 'hmmtop2')) {
					$this_method_topology =~ s/H/M/g;
				} elsif ($method eq 'phdthtm_pbil') {
					$this_method_topology =~ s/T/M/g;
				} elsif ($method eq 'toppred2') {
					$this_method_topology =~ s/l/\_/g;
				} elsif ($method eq 'memsat3') {
					$this_method_topology =~ s/H/M/g;
					$this_method_topology =~ s/\+/i/g;
					$this_method_topology =~ s/\-/o/g;
					$this_method_topology_with_interface_as_nonmembrane = $this_method_topology;
					$this_method_topology =~ s/I/M/g;
					$this_method_topology =~ s/O/M/g;
					$this_method_topology_with_interface_as_nonmembrane =~ s/I/i/g;
					$this_method_topology_with_interface_as_nonmembrane =~ s/O/o/g;
				} elsif ($method eq 'octopus') {
					$this_method_topology =~ s/L/M/g;
					$this_method_topology =~ s/R/M/g;
					$this_method_topology =~ s/g/\?/g; # is an extension of whatever inside/outside precedes/follows it
					$this_method_topology_with_interface_as_nonmembrane = $this_method_topology;
					$this_method_topology =~ s/I/M/g;
					$this_method_topology_with_interface_as_nonmembrane =~ s/I/\?/g; # is an extension of whatever inside/outside precedes/follows it
					if ( index($this_method_topology, '?') > -1 ) {
						my @this_topology_array = split( //, $this_method_topology );
						$this_method_topology = '';
						my $prev_nonmembrane_char = '_';
						for ( my $j = length(@this_topology_array) - 1; $j >= 0; $j-- ) {
							if ($this_topology_array[$j] eq 'M') {
								$prev_nonmembrane_char = '_';
							} elsif (($this_topology_array[$j] eq 'i') or ($this_topology_array[$j] eq 'o')) {
								$prev_nonmembrane_char = $this_topology_array[$j];
							} elsif (($this_topology_array[$j] eq '?') and ($prev_nonmembrane_char ne '_')) {
								$this_topology_array[$j] = $prev_nonmembrane_char;
							}
						}
						for ( my $j = 0; $j < length(@this_topology_array); $j++ ) {
							if ($this_topology_array[$j] eq 'M') {
								$prev_nonmembrane_char = '_';
							} elsif (($this_topology_array[$j] eq 'i') or ($this_topology_array[$j] eq 'o')) {
								$prev_nonmembrane_char = $this_topology_array[$j];
							} elsif (($this_topology_array[$j] eq '?') and ($prev_nonmembrane_char ne '_')) {
								$this_topology_array[$j] = $prev_nonmembrane_char;
							}
							$this_method_topology .= $this_topology_array[$j];
						}
					}
					if (($method eq 'octopus') or ($method eq 'memsat3')) {
						if ( index($this_method_topology_with_interface_as_nonmembrane, '?') > -1 ) {
							my @this_topology_array = split( //, $this_method_topology_with_interface_as_nonmembrane );
							$this_method_topology_with_interface_as_nonmembrane = '';
							my $prev_nonmembrane_char = '_';
							for ( my $j = length(@this_topology_array) - 1; $j >= 0; $j-- ) {
								if ($this_topology_array[$j] eq 'M') {
									$prev_nonmembrane_char = '_';
								} elsif (($this_topology_array[$j] eq 'i') or ($this_topology_array[$j] eq 'o')) {
									$prev_nonmembrane_char = $this_topology_array[$j];
								} elsif (($this_topology_array[$j] eq '?') and ($prev_nonmembrane_char ne '_')) {
									$this_topology_array[$j] = $prev_nonmembrane_char;
								}
							}
							for ( my $j = 0; $j < length(@this_topology_array); $j++ ) {
								if ($this_topology_array[$j] eq 'M') {
									$prev_nonmembrane_char = '_';
								} elsif (($this_topology_array[$j] eq 'i') or ($this_topology_array[$j] eq 'o')) {
									$prev_nonmembrane_char = $this_topology_array[$j];
								} elsif (($this_topology_array[$j] eq '?') and ($prev_nonmembrane_char ne '_')) {
									$this_topology_array[$j] = $prev_nonmembrane_char;
								}
								$this_method_topology_with_interface_as_nonmembrane .= $this_topology_array[$j];
							}
						}
					}
				}

				my @shifted_predicted_topology_residues = split( //, $this_method_topology);
				for ( my $i = 1; $i <= $stats_num_residues; $i++ ) {
					$predicted_topology_residues[$i] = $shifted_predicted_topology_residues[ $i - 1 ];
				}
				my $prev_char = '';
				my $prev_start = 0;
				my $prev_end = 0;
				for ( my $i = 1; $i <= $stats_num_residues; $i++ ) {
					my $this_char = $predicted_topology_residues[$i];
					if ($this_char eq '_') {
						$this_char = $prev_char;
					}
					if ($this_char ne $prev_char) {
						# process the previous segment
						if ( ($prev_start != 0) and (($prev_char eq 'i') or ($prev_char eq 'o')) ) {
							$prev_end = $i - 1;
							push( @predicted_topology_start, $prev_start );
							push( @predicted_topology_end, $prev_end );
							push( @predicted_topology_location, $prev_char );
						}
						# start processing of this new segment
						if (($this_char eq 'i') or ($this_char eq 'o')) {
							$prev_start = $i;
						}
					}
					$prev_char = $this_char;
				}
				# process the last segment
				if ( ($prev_start != 0) and (($prev_char eq 'i') or ($prev_char eq 'o')) ) {
					$prev_end = $stats_num_residues;
					push( @predicted_topology_start, $prev_start );
					push( @predicted_topology_end, $prev_end );
					push( @predicted_topology_location, $prev_char );

				}
				$stats_num_prd_iosegment_in_this_seq = @predicted_topology_start;
				$stats_num_prd_iosegment_in_this_method += $stats_num_prd_iosegment_in_this_seq;

				if (($method eq 'octopus') or ($method eq 'memsat3')) {
					my @shifted_predicted_topology_residues_with_interface_as_nonmembrane = split( //, $this_method_topology_with_interface_as_nonmembrane);
					for ( my $i = 1; $i <= $stats_num_residues; $i++ ) {
						$predicted_topology_residues_with_interface_as_nonmembrane[$i] = $shifted_predicted_topology_residues_with_interface_as_nonmembrane[ $i - 1 ];
					}
					$prev_char = '';
					$prev_start = 0;
					$prev_end = 0;
					for ( my $i = 1; $i <= $stats_num_residues; $i++ ) {
						my $this_char = $predicted_topology_residues_with_interface_as_nonmembrane[$i];
						if ($this_char eq '_') {
							$this_char = $prev_char;
						}
						if ($this_char ne $prev_char) {
							# process the previous segment
							if ( ($prev_start != 0) and (($prev_char eq 'i') or ($prev_char eq 'o')) ) {
								$prev_end = $i - 1;
								push( @predicted_topology_start_with_interface_as_nonmembrane, $prev_start );
								push( @predicted_topology_end_with_interface_as_nonmembrane, $prev_end );
								push( @predicted_topology_location_with_interface_as_nonmembrane, $prev_char );
							}
							# start processing of this new segment
							if (($this_char eq 'i') or ($this_char eq 'o')) {
								$prev_start = $i;
							}
						}
						$prev_char = $this_char;
					}
					# process the last segment
					if ( ($prev_start != 0) and (($prev_char eq 'i') or ($prev_char eq 'o')) ) {
						$prev_end = $stats_num_residues;
						push( @predicted_topology_start_with_interface_as_nonmembrane, $prev_start );
						push( @predicted_topology_end_with_interface_as_nonmembrane, $prev_end );
						push( @predicted_topology_location_with_interface_as_nonmembrane, $prev_char );
					}
				}
			}
		}

		# for each observed membrane segment in this sequence, find out if it was predicted
		for ( my $i = 0; $i < @observed_helices; $i++ ) {
			my $observed_start = $observed_helices_start[$i];
			my $observed_end = $observed_helices_end[$i];
			my $found_overlapping_tms = 0;
			my $found_a_predicted_tms_for_this_observed_tms = 0;
			my $minimum_overlap_for_this_segment = $global->{'minimum_overlap'};
			my $length_of_this_segment = $observed_end - $observed_start + 1;
			if ($length_of_this_segment < $minimum_overlap_for_this_segment) {
				$minimum_overlap_for_this_segment = $length_of_this_segment;
			}
			for ( my $j = 0; $j < @predicted_helices; $j++ ) {
				my $predicted_start = $predicted_helices_start[$j];
				my $predicted_end = $predicted_helices_end[$j];
				if (($found_overlapping_tms == 0) && ($predicted_helices_already_counted[$j] == 0)) {
					if ( (($observed_start <= $predicted_start) && ($predicted_start <= $observed_end)) ||
						(($observed_start <= $predicted_end) && ($predicted_end <= $observed_end)) ||
						(($predicted_start <= $observed_start) && ($observed_end <= $predicted_end)) ) {
						my $num_overlap = 0;
						for ( my $k = $observed_start; $k <= $observed_end; $k++ ) {
							if (($k >= 1) && ($k <= $stats_num_residues)) {
								if ($predicted_residues[$k] == 1) {
									$num_overlap++;
								}
							}
						}
						if ($num_overlap >= $minimum_overlap_for_this_segment) {
							$stats_num_correct_obs_tms_in_this_seq++;
							$stats_num_correct_obs_tms_in_this_method = $stats_num_correct_obs_tms_in_this_method + 1;
							$found_a_predicted_tms_for_this_observed_tms = 1;
							$predicted_helices_already_counted[$j] = 1;
							$found_overlapping_tms = 1;
							# calculate the sum to use for calculation of average helix end difference in residues
							my $this_diff_start = abs($observed_start - $predicted_start);
							my $this_diff_end = abs($observed_end - $predicted_end);
							$stats_sum_dist_obs_helix_end_minus_prd_helix_end = $stats_sum_dist_obs_helix_end_minus_prd_helix_end + $this_diff_start + $this_diff_end;
						}

						my $start_dist = abs($observed_start - $predicted_start);
						my $end_dist = abs($observed_end - $predicted_end);

						if ($start_dist <= $global->{'max_dist_approx_ends'}) {
							$stats_num_correct_obs_approx_ends_in_this_method = $stats_num_correct_obs_approx_ends_in_this_method + 1;
						}
						if ($end_dist <= $global->{'max_dist_approx_ends'}) {
							$stats_num_correct_obs_approx_ends_in_this_method = $stats_num_correct_obs_approx_ends_in_this_method + 1;
						}
					}
				}
			}
		}

		# for each observed topology segment in this sequence, find out if it was predicted
		if (($this_method_has_topology == 1) and ($this_protein_has_topology == 1)) {
			my $seen_nterminal_topology = 0;
			for ( my $i = 0; $i < @observed_topology_start; $i++ ) {
				my $observed_start = $observed_topology_start[$i];
				my $observed_end = $observed_topology_end[$i];
				my $num_predicted_i = 0;
				my $num_predicted_o = 0;
				for ( my $k = $observed_start; $k <= $observed_end; $k++ ) {
					if (($predicted_topology_residues[$k] eq 'i') or ($predicted_topology_residues[$k] eq 'I') or ($predicted_topology_residues[$k] eq '+')) {
						$num_predicted_i++;
					} elsif (($predicted_topology_residues[$k] eq 'o') or ($predicted_topology_residues[$k] eq 'O') or ($predicted_topology_residues[$k] eq '-')) {
						$num_predicted_o++;
					}
				}
				if (($num_predicted_i >= $num_predicted_o) && ($num_predicted_i > 0)) {
					if ($observed_topology_location[$i] eq 'i') {
						$stats_num_correct_obs_iosegment_in_this_seq += 1;
						$stats_num_correct_obs_iosegment_in_this_method += 1;
						if ($seen_nterminal_topology == 0) {
							$stats_num_correct_nterminal_topology += 1;
						}
					}
				} elsif (($num_predicted_o >= $num_predicted_i) && ($num_predicted_o > 0)) {
					if ($observed_topology_location[$i] eq 'o') {
						$stats_num_correct_obs_iosegment_in_this_seq += 1;
						$stats_num_correct_obs_iosegment_in_this_method += 1;
						if ($seen_nterminal_topology == 0) {
							$stats_num_correct_nterminal_topology += 1;
						}
					}
				}
				# if the observed i/o segment didn't quite overlap with a predicted i/o segment, look a little further on either side
				if (($num_predicted_i == 0) && ($num_predicted_o == 0)) {
					for ( my $k = ($observed_start - 5); $k <= ($observed_end + 5); $k++ ) {
						if (($k >= 0) && ($k < $stats_num_residues)) {
							if (($predicted_topology_residues[$k] eq 'i') or ($predicted_topology_residues[$k] eq 'I') or ($predicted_topology_residues[$k] eq '+')) {
								$num_predicted_i++;
							} elsif (($predicted_topology_residues[$k] eq 'o') or ($predicted_topology_residues[$k] eq 'O') or ($predicted_topology_residues[$k] eq '-')) {
								$num_predicted_o++;
							}
						}
					}
					if (($num_predicted_i >= $num_predicted_o) && ($num_predicted_i > 0)) {
						if ($observed_topology_location[$i] eq 'i') {
							$stats_num_correct_obs_iosegment_in_this_seq += 1;
							$stats_num_correct_obs_iosegment_in_this_method += 1;
							if ($seen_nterminal_topology == 0) {
								$stats_num_correct_nterminal_topology += 1;
							}
						}
					} elsif (($num_predicted_o >= $num_predicted_i) && ($num_predicted_o > 0)) {
						if ($observed_topology_location[$i] eq 'o') {
							$stats_num_correct_obs_iosegment_in_this_seq += 1;
							$stats_num_correct_obs_iosegment_in_this_method += 1;
							if ($seen_nterminal_topology == 0) {
								$stats_num_correct_nterminal_topology += 1;
							}
						}
					}
				}
				$seen_nterminal_topology = 1;
			}
			for ( my $i = 0; $i < @predicted_topology_start; $i++ ) {
				my $predicted_start = $predicted_topology_start[$i];
				my $predicted_end = $predicted_topology_end[$i];
				for ( my $k = $predicted_start; $k <= $predicted_end; $k++ ) {
					if ($predicted_topology_residues[$k] eq $observed_topology_residues[$k]) {
						if ($predicted_topology_residues[$k] eq 'i') {
							$stats_num_correct_prd_i_topology_residues_in_this_method++;
						} elsif ($predicted_topology_residues[$k] eq 'o') {
							$stats_num_correct_prd_o_topology_residues_in_this_method++;
						}
					}
				}
			}
		}

		# for each predicted membrane segment in this sequence, find out if it was observed
		for ( my $i = 0; $i < @predicted_helices; $i++ ) {
			my $predicted_start = $predicted_helices_start[$i];
			my $predicted_end = $predicted_helices_end[$i];
			my $found_overlapping_tms = 0;
			my $found_an_observed_tms_for_this_predicted_tms = 0;
			my $minimum_overlap_for_this_segment = $global->{'minimum_overlap'};
			my $length_of_this_segment = $predicted_end - $predicted_start + 1;
			if ($length_of_this_segment < $minimum_overlap_for_this_segment) {
				$minimum_overlap_for_this_segment = $length_of_this_segment;
			}
			for ( my $j = 0; $j < @observed_helices; $j++ ) {
				my $observed_start = $observed_helices_start[$j];
				my $observed_end = $observed_helices_end[$j];
				if (($found_overlapping_tms == 0) && ($observed_helices_already_counted[$j] == 0)) {
					if ( (($predicted_start <= $observed_start) && ($observed_start <= $predicted_end)) ||
						(($predicted_start <= $observed_end) && ($observed_end <= $predicted_end)) ||
						(($observed_start <= $predicted_start) && ($predicted_end <= $observed_end)) ) {
						my $num_overlap = 0;
						for ( my $k = $predicted_start; $k <= $predicted_end; $k++ ) {
							if (($k >= 1) && ($k <= $stats_num_residues)) {
								if ($observed_residues[$k] == 1) {
									$num_overlap++;
								}
							}
						}
						if ($num_overlap >= $minimum_overlap_for_this_segment) {
							$stats_num_correct_prd_tms_in_this_seq++;
							$stats_num_correct_prd_tms_in_this_method = $stats_num_correct_prd_tms_in_this_method + 1;
							$found_an_observed_tms_for_this_predicted_tms = 1;
							$observed_helices_already_counted[$j] = 1;
							$found_overlapping_tms = 1;
						}

						my $start_dist = abs($observed_start - $predicted_start);
						my $end_dist = abs($observed_end - $predicted_end);

						if ($start_dist <= $global->{'max_dist_approx_ends'}) {
							$stats_num_correct_prd_approx_ends_in_this_method = $stats_num_correct_prd_approx_ends_in_this_method + 1;
						}
						if ($end_dist <= $global->{'max_dist_approx_ends'}) {
							$stats_num_correct_prd_approx_ends_in_this_method = $stats_num_correct_prd_approx_ends_in_this_method + 1;
						}

						if ($start_dist <= $global->{'variance'}) {
							my $distance = $start_dist;
							if ($distance == 0) {
								$stats_sum_of_scaled_scores_of_helix_ends = $stats_sum_of_scaled_scores_of_helix_ends + 1;
							} else {
								my $distance_squared = $distance * $distance;
								my $score_for_this_end_1 = exp( -1 * ($distance_squared / $global->{'variance2'}) );
								my $score_for_this_end_2 = $score_for_this_end_1 / $global->{'std_dev_root_2_pi'};
								my $score_for_this_end_3 = $score_for_this_end_2 / $global->{'helix_ends_score_denominator'};
								$stats_sum_of_scaled_scores_of_helix_ends = $stats_sum_of_scaled_scores_of_helix_ends + $score_for_this_end_3;
							}
						}
						if ($end_dist <= $global->{'variance'}) {
							my $distance = $end_dist;
							if ($distance == 0) {
								$stats_sum_of_scaled_scores_of_helix_ends = $stats_sum_of_scaled_scores_of_helix_ends + 1;
							} else {
								my $distance_squared = $distance * $distance;
								my $score_for_this_end_1 = exp( -1 * ($distance_squared / $global->{'variance2'}) );
								my $score_for_this_end_2 = $score_for_this_end_1 / $global->{'std_dev_root_2_pi'};
								my $score_for_this_end_3 = $score_for_this_end_2 / $global->{'helix_ends_score_denominator'};
								$stats_sum_of_scaled_scores_of_helix_ends = $stats_sum_of_scaled_scores_of_helix_ends + $score_for_this_end_3;
							}
						}
					}
				}
			}
		}

		# for each predicted topology segment in this sequence, find out if it was observed
		if (($this_method_has_topology == 1) and ($this_protein_has_topology == 1)) {
			if (($method eq 'octopus') or ($method eq 'memsat3')) {
				for ( my $i = 0; $i < @predicted_topology_start_with_interface_as_nonmembrane; $i++ ) {
					my $predicted_start = $predicted_topology_start_with_interface_as_nonmembrane[$i];
					my $predicted_end = $predicted_topology_end_with_interface_as_nonmembrane[$i];
					my $found_overlapping_segment = 0;
					my $found_an_observed_segment_for_this_predicted_segment = 0;
					my $minimum_overlap_for_this_segment = $global->{'minimum_overlap'};
					my $length_of_this_segment = $predicted_end - $predicted_start + 1;
					if ($length_of_this_segment < $minimum_overlap_for_this_segment) {
						$minimum_overlap_for_this_segment = $length_of_this_segment;
					}
					for ( my $j = 0; $j < @observed_topology_start; $j++ ) {
						my $observed_start = $observed_topology_start[$j];
						my $observed_end = $observed_topology_end[$j];
						if (($found_overlapping_segment == 0) && ($observed_topology_already_counted[$j] == 0)) {
							if ( (($predicted_start <= $observed_start) && ($observed_start <= $predicted_end)) ||
								(($predicted_start <= $observed_end) && ($observed_end <= $predicted_end)) ||
								(($observed_start <= $predicted_start) && ($predicted_end <= $observed_end)) ) {
								my $num_overlap = 0;
								for ( my $k = $predicted_start; $k <= $predicted_end; $k++ ) {
									# if this observed segment is a membrane segment, then prediction interface residues will be counted as membrane residues.
									# if this observed segment is a non-membrane segment (i or o), 
									# then prediction interface residues will be counted as whatever non-membrane (i or o) residues that they are next to.
									my $predicted_topology_location_to_use_for_compare = $predicted_topology_location_with_interface_as_nonmembrane[$i];
									if ($observed_topology_location[$j] eq 'M') {
										$predicted_topology_location_to_use_for_compare = $predicted_topology_location[$i];
									}
									if (($k >= 1) && ($k <= $stats_num_residues)) {
										if ($observed_topology_location[$j] eq $predicted_topology_location_to_use_for_compare) {
											$num_overlap++;
										}
									}
								}
								if ($num_overlap >= $minimum_overlap_for_this_segment) {
									$stats_num_correct_prd_iosegment_in_this_seq++;
									$stats_num_correct_prd_iosegment_in_this_method = $stats_num_correct_prd_iosegment_in_this_method + 1;
									$found_an_observed_segment_for_this_predicted_segment = 1;
									$observed_topology_already_counted[$j] = 1;
									$found_overlapping_segment = 1;
								}
							}
						}
					}
				}
			} else { # the rest of the methods don't have this complication caused by predicted interface residues
				for ( my $i = 0; $i < @predicted_topology_start; $i++ ) {
					my $predicted_start = $predicted_topology_start[$i];
					my $predicted_end = $predicted_topology_end[$i];
					my $found_overlapping_segment = 0;
					my $found_an_observed_segment_for_this_predicted_segment = 0;
					my $minimum_overlap_for_this_segment = $global->{'minimum_overlap'};
					my $length_of_this_segment = $predicted_end - $predicted_start + 1;
					if ($length_of_this_segment < $minimum_overlap_for_this_segment) {
						$minimum_overlap_for_this_segment = $length_of_this_segment;
					}
					for ( my $j = 0; $j < @observed_topology_start; $j++ ) {
						my $observed_start = $observed_topology_start[$j];
						my $observed_end = $observed_topology_end[$j];
						if (($found_overlapping_segment == 0) && ($observed_topology_already_counted[$j] == 0)) {
							if ( (($predicted_start <= $observed_start) && ($observed_start <= $predicted_end)) ||
								(($predicted_start <= $observed_end) && ($observed_end <= $predicted_end)) ||
								(($observed_start <= $predicted_start) && ($predicted_end <= $observed_end)) ) {
								my $num_overlap = 0;
								for ( my $k = $predicted_start; $k <= $predicted_end; $k++ ) {
									if (($k >= 1) && ($k <= $stats_num_residues)) {
										if ($observed_topology_location[$j] eq $predicted_topology_location[$i]) {
											$num_overlap++;
										}
									}
								}
								if ($num_overlap >= $minimum_overlap_for_this_segment) {
									$stats_num_correct_prd_iosegment_in_this_seq++;
									$stats_num_correct_prd_iosegment_in_this_method = $stats_num_correct_prd_iosegment_in_this_method + 1;
									$found_an_observed_segment_for_this_predicted_segment = 1;
									$observed_topology_already_counted[$j] = 1;
									$found_overlapping_segment = 1;
								}
							}
						}
					}
				}
			}
		}

		# find out if the segments predicted for this sequence matches the observed segments
		if (($stats_num_obs_tms_in_this_seq == $stats_num_prd_tms_in_this_seq) && ($stats_num_obs_tms_in_this_seq == $stats_num_correct_prd_tms_in_this_seq)) {
			$stats_qok_for_this_seq = 1;
			$stats_num_correct_obs_and_prd_tms_in_this_method += 1;
		}

		# find out if the topology prediction for this sequence matches the observed topoogy
		if (($this_method_has_topology == 1) and ($this_protein_has_topology == 1)) {
			if ($stats_num_correct_obs_iosegment_in_this_seq == $stats_num_obs_iosegment_in_this_seq) {
				$stats_qok3_for_this_seq = 1;
				$stats_qok3_count_for_this_method += 1;
				$stats_num_correct_obs_and_prd_iosegment_in_this_method += 1;
			}
			if ($stats_num_obs_iosegment_in_this_seq > 0) {
				$stats_ratio_correct_iotopology_segments_in_this_method += $stats_num_correct_obs_iosegment_in_this_seq / $stats_num_obs_iosegment_in_this_seq;
			}
		}

		# count how many residues are correctly predicted as being tms or not tms

		my $stats_num_obs_tms_residues_in_this_seq = 0;
		my $stats_num_prd_tms_residues_in_this_seq = 0;
		my $stats_num_correct_obs_tms_residues_in_this_seq = 0;
		my $stats_num_correct_prd_tms_residues_in_this_seq = 0;
		my $stats_num_obs_nontms_residues_in_this_seq = 0;
		my $stats_num_prd_nontms_residues_in_this_seq = 0;
		my $stats_num_correct_obs_nontms_residues_in_this_seq = 0;
		my $stats_num_correct_prd_nontms_residues_in_this_seq = 0;
		my $stats_num_correct_residues_in_this_seq = 0;

		for ( my $i = 1; $i <= $stats_num_residues; $i++ ) {
			if ($observed_residues[$i] == 1) {
				$stats_num_obs_tms_residues_in_this_seq++;
				if ($predicted_residues[$i] == 1) {
					$stats_num_correct_obs_tms_residues_in_this_seq++;
				}
			} else { # $observed_residues[$i] == 0
				$stats_num_obs_nontms_residues_in_this_seq++;
				if ($predicted_residues[$i] == 0) {
					$stats_num_correct_obs_nontms_residues_in_this_seq++;
				}
			}
			if ($predicted_residues[$i] == 1) {
				$stats_num_prd_tms_residues_in_this_seq++;
				if ($observed_residues[$i] == 1) {
					$stats_num_correct_prd_tms_residues_in_this_seq++;
				}
			} else { # $predicted_residues[$i] == 0
				$stats_num_prd_nontms_residues_in_this_seq++;
				if ($observed_residues[$i] == 0) {
					$stats_num_correct_prd_nontms_residues_in_this_seq++;
				}
			}
			if (($observed_residues[$i] == 1) && ($predicted_residues[$i] == 1)) {
				$stats_num_correct_residues_in_this_seq++;
			} elsif (($observed_residues[$i] == 0) && ($predicted_residues[$i] == 0)) {
				$stats_num_correct_residues_in_this_seq++;
			}
		}

		# add the counts of correct residues to the overall statistics

		if ($stats_num_residues > 0) {
			$stats_ratio_correct_residues_in_this_method += ($stats_num_correct_residues_in_this_seq / $stats_num_residues);
		}
		if ($stats_num_obs_tms_residues_in_this_seq > 0) {
			$stats_num_correct_obs_tms_residues_in_this_method += $stats_num_correct_obs_tms_residues_in_this_seq;
			$stats_num_obs_tms_residues_in_this_method += $stats_num_obs_tms_residues_in_this_seq;
		}
		if ($stats_num_prd_tms_residues_in_this_seq > 0) {
			$stats_num_correct_prd_tms_residues_in_this_method += $stats_num_correct_prd_tms_residues_in_this_seq;
			$stats_num_prd_tms_residues_in_this_method += $stats_num_prd_tms_residues_in_this_seq;
		}
		if ($stats_num_obs_nontms_residues_in_this_seq > 0) {
			$stats_num_correct_obs_nontms_residues_in_this_method += $stats_num_correct_obs_nontms_residues_in_this_seq;
			$stats_num_obs_nontms_residues_in_this_method += $stats_num_obs_nontms_residues_in_this_seq;
		}
		if ($stats_num_prd_nontms_residues_in_this_seq > 0) {
			$stats_num_correct_prd_nontms_residues_in_this_method += $stats_num_correct_prd_nontms_residues_in_this_seq;
			$stats_num_prd_nontms_residues_in_this_method += $stats_num_prd_nontms_residues_in_this_seq;
		}

		# count how many residues have correctly predicted topology

		my $stats_num_correct_iotopology_residues_in_this_seq = 0;
		my $stats_num_correct_iomtopology_residues_in_this_seq = 0;

		if (($this_method_has_topology == 1) and ($this_protein_has_topology == 1)) {
			if (($method eq 'octopus') or ($method eq 'memsat3')) {
				for ( my $i = 1; $i <= $stats_num_residues; $i++ ) {
					# if this observed segment is a membrane segment, then prediction interface residues will be counted as membrane residues.
					# if this observed segment is a non-membrane segment (i or o), 
					# then prediction interface residues will be counted as whatever non-membrane (i or o) residues that they are next to.
					my $predicted_topology_residues_to_use_for_compare = $predicted_topology_residues_with_interface_as_nonmembrane[$i];
					if ($observed_topology_residues[$i] eq 'M') {
						$predicted_topology_residues_to_use_for_compare = $predicted_topology_residues[$i];
					}
					if ($observed_topology_residues[$i] eq $predicted_topology_residues_to_use_for_compare) {
						if (($observed_topology_residues[$i] eq 'i') or ($observed_topology_residues[$i] eq 'o')) {
							$stats_num_correct_iotopology_residues_in_this_seq++;
						}
						if (($observed_topology_residues[$i] eq 'i') or ($observed_topology_residues[$i] eq 'o') or ($observed_topology_residues[$i] eq 'M')) {
							$stats_num_correct_iomtopology_residues_in_this_seq++;
						}
					}
				}
			} else { # the rest of the methods don't have this complication caused by predicted interface residues
				for ( my $i = 1; $i <= $stats_num_residues; $i++ ) {
					if ($observed_topology_residues[$i] eq $predicted_topology_residues[$i]) {
						if (($observed_topology_residues[$i] eq 'i') or ($observed_topology_residues[$i] eq 'o')) {
							$stats_num_correct_iotopology_residues_in_this_seq++;
						}
						if (($observed_topology_residues[$i] eq 'i') or ($observed_topology_residues[$i] eq 'o') or ($observed_topology_residues[$i] eq 'M')) {
							$stats_num_correct_iomtopology_residues_in_this_seq++;
						}
					}
				}
			}
		}

		# add the counts of correct topology residues to the overall statistics

		if ($stats_num_known_topology_residues > 0) {
			$stats_ratio_correct_iotopology_residues_in_this_method += ($stats_num_correct_iotopology_residues_in_this_seq / $stats_num_known_io_topology_residues);
			$stats_ratio_correct_iomtopology_residues_in_this_method += ($stats_num_correct_iomtopology_residues_in_this_seq / $stats_num_known_topology_residues);
		}
	}

	# do the final calculations that will be displayed for this method's statistics, and format the figures for display

	if ($stats_num_seqs > 0 ) {
		$stats->{'qok'} = $stats_num_correct_obs_and_prd_tms_in_this_method * 100 / $stats_num_seqs;
		$stats->{'qok'} = sprintf("%0.f", $stats->{'qok'});
	} else {
		$stats->{'qok'} = 0;
	}

	if ($stats_num_obs_tms_in_this_method > 0) {
		$stats->{'qhtm_obs'} = $stats_num_correct_obs_tms_in_this_method * 100 / $stats_num_obs_tms_in_this_method;
		$stats->{'qhtm_obs'} = sprintf("%0.f", $stats->{'qhtm_obs'});
	} else {
		if ($stats_num_correct_obs_tms_in_this_method == 0) {
			$stats->{'qhtm_obs'} = 0;
		} else {
			$stats->{'qhtm_obs'} = 100;
		}
	}

	if ($stats_num_prd_tms_in_this_method > 0) {
		$stats->{'qhtm_prd'} = $stats_num_correct_prd_tms_in_this_method * 100 / $stats_num_prd_tms_in_this_method;
		$stats->{'qhtm_prd'} = sprintf("%0.f", $stats->{'qhtm_prd'});
	} else {
		if ($stats_num_correct_prd_tms_in_this_method == 0) {
			$stats->{'qhtm_prd'} = 0;
		} else {
			$stats->{'qhtm_prd'} = 100;
		}
	}

	if ($stats_num_seqs > 0 ) {
		$stats->{'q2'} = $stats_ratio_correct_residues_in_this_method * 100 / $stats_num_seqs;
		$stats->{'q2'} = sprintf("%0.f", $stats->{'q2'});
	} else {
		$stats->{'q2'} = 0;
	}

	if ($stats_num_obs_tms_residues_in_this_method > 0) {
		$stats->{'q2t_obs'} = ($stats_num_correct_prd_tms_residues_in_this_method * 100) / $stats_num_obs_tms_residues_in_this_method;
	} else {
		if ($stats_num_correct_prd_tms_residues_in_this_method == 0) {
			$stats->{'q2t_obs'} = 0;
		} else {
			$stats->{'q2t_obs'} = 100;
		}
	}
	$stats->{'q2t_obs'} = sprintf("%0.f", $stats->{'q2t_obs'});

	if ($stats_num_prd_tms_residues_in_this_method > 0) {
		$stats->{'q2t_prd'} = ($stats_num_correct_prd_tms_residues_in_this_method * 100) / $stats_num_prd_tms_residues_in_this_method;
	} else {
		if ($stats_num_correct_prd_tms_residues_in_this_method == 0) {
			$stats->{'q2t_prd'} = 0;
		} else {
			$stats->{'q2t_prd'} = 100;
		}
	}
	$stats->{'q2t_prd'} = sprintf("%0.f", $stats->{'q2t_prd'});

	if ($stats_num_obs_nontms_residues_in_this_method > 0) {
		$stats->{'q2n_obs'} = ($stats_num_correct_prd_nontms_residues_in_this_method * 100) / $stats_num_obs_nontms_residues_in_this_method;
	} else {
		if ($stats_num_correct_prd_nontms_residues_in_this_method == 0) {
			$stats->{'q2n_obs'} = 0;
		} else {
			$stats->{'q2n_obs'} = 100;
		}
	}
	$stats->{'q2n_obs'} = sprintf("%0.f", $stats->{'q2n_obs'});

	if ($stats_num_prd_nontms_residues_in_this_method > 0) {
		$stats->{'q2n_prd'} = ($stats_num_correct_prd_nontms_residues_in_this_method * 100) / $stats_num_prd_nontms_residues_in_this_method;
	} else {
		if ($stats_num_correct_prd_nontms_residues_in_this_method == 0) {
			$stats->{'q2n_prd'} = 0;
		} else {
			$stats->{'q2n_prd'} = 100;
		}
	}
	$stats->{'q2n_prd'} = sprintf("%0.f", $stats->{'q2n_prd'});

	if ($stats_num_obs_approx_ends_in_this_method > 0) {
		$stats->{'qhe_obs'} = $stats_num_correct_obs_approx_ends_in_this_method * 100 / $stats_num_obs_approx_ends_in_this_method;
		$stats->{'qhe_obs'} = sprintf("%0.f", $stats->{'qhe_obs'});
	} else {
		if ($stats_num_correct_obs_approx_ends_in_this_method == 0) {
			$stats->{'qhe_obs'} = 0;
		} else {
			$stats->{'qhe_obs'} = 100;
		}
	}

	# calculate the average helix end difference
	if ($stats_num_correct_obs_tms_in_this_method > 0) {
		$stats->{'avhe_diff'} = $stats_sum_dist_obs_helix_end_minus_prd_helix_end / ($stats_num_correct_obs_tms_in_this_method * 2);
	} else {
		$stats->{'avhe_diff'} = 0;
	}

	if ($stats_num_prd_approx_ends_in_this_method > 0) {
		$stats->{'qhe_prd'} = $stats_num_correct_prd_approx_ends_in_this_method * 100 / $stats_num_prd_approx_ends_in_this_method;
		$stats->{'qhe_prd'} = sprintf("%0.f", $stats->{'qhe_prd'});
	} else {
		if ($stats_num_correct_prd_approx_ends_in_this_method == 0) {
			$stats->{'qhe_prd'} = 0;
		} else {
			$stats->{'qhe_prd'} = 100;
		}
	}

	if ($stats_num_obs_approx_ends_in_this_method > 0) {
		$stats->{'qnhe_obs'} = $stats_sum_of_scaled_scores_of_helix_ends * 100 / $stats_num_obs_approx_ends_in_this_method;
		$stats->{'qnhe_obs'} = sprintf("%0.f", $stats->{'qnhe_obs'});
	} else {
		if ($stats_num_obs_approx_ends_in_this_method == 0) {
			$stats->{'qnhe_obs'} = 0;
		} else {
			$stats->{'qnhe_obs'} = 100;
		}
	}

	if ($stats_num_prd_approx_ends_in_this_method > 0) {
		$stats->{'qnhe_prd'} = $stats_sum_of_scaled_scores_of_helix_ends * 100 / $stats_num_prd_approx_ends_in_this_method;
		$stats->{'qnhe_prd'} = sprintf("%0.f", $stats->{'qnhe_prd'});
	} else {
		if ($stats_sum_of_scaled_scores_of_helix_ends == 0) {
			$stats->{'qnhe_prd'} = 0;
		} else {
			$stats->{'qnhe_prd'} = 100;
		}
	}

	my $htm_mcc_TP = $stats_num_correct_prd_tms_residues_in_this_method;
	my $htm_mcc_TN = $stats_num_correct_prd_nontms_residues_in_this_method;
	my $htm_mcc_FP = $stats_num_obs_tms_residues_in_this_method - $stats_num_correct_prd_tms_residues_in_this_method;
	my $htm_mcc_FN = $stats_num_obs_nontms_residues_in_this_method - $stats_num_correct_prd_nontms_residues_in_this_method;
	my $htm_mcc_top = ($htm_mcc_TP * $htm_mcc_TN) - ($htm_mcc_FP * $htm_mcc_FN);
	my $htm_mcc_bottom = sqrt( ($htm_mcc_TP + $htm_mcc_FP) * ($htm_mcc_TP + $htm_mcc_FN) * ($htm_mcc_TN + $htm_mcc_FP) * ($htm_mcc_TN + $htm_mcc_FN) );
	if ($htm_mcc_bottom != 0) {
		$stats->{'htm_mcc'} = $htm_mcc_top / $htm_mcc_bottom;
	}

	if ($this_method_has_topology == 1) {

		if ($stats_num_seqs > 0 ) {
			$stats->{'qok3'} = $stats_qok3_count_for_this_method * 100 / $stats_num_seqs;
			$stats->{'qok3'} = sprintf("%0.f", $stats->{'qok3'});
		} else {
			$stats->{'qok3'} = 0;
		}

		if ($stats_num_seqs > 0 ) {
			$stats->{'nterm'} = $stats_num_correct_nterminal_topology * 100 / $stats_num_seqs;
		} else {
			$stats->{'nterm'} = 0;
		}

		if ($stats_num_seqs > 0 ) {
			$stats->{'q2_ioseg'} = $stats_ratio_correct_iotopology_segments_in_this_method * 100 / $stats_num_seqs;
			$stats->{'q2_ioseg'} = sprintf("%0.f", $stats->{'q2_ioseg'});
		} else {
			$stats->{'q2_ioseg'} = 0;
		}

		$stats_num_obs_iomsegment_in_this_method = $stats_num_obs_tms_in_this_method + $stats_num_obs_iosegment_in_this_method;
		$stats_num_correct_obs_iomsegment_in_this_method = $stats_num_correct_obs_tms_in_this_method + $stats_num_correct_obs_iosegment_in_this_method;
		if ($stats_num_obs_iomsegment_in_this_method > 0) {
			$stats->{'qiom_obs'} = $stats_num_correct_obs_iomsegment_in_this_method * 100 / $stats_num_obs_iomsegment_in_this_method;
			$stats->{'qiom_obs'} = sprintf("%0.f", $stats->{'qiom_obs'});
		} else {
			if ($stats_num_correct_obs_iomsegment_in_this_method == 0) {
				$stats->{'qiom_obs'} = 0;
			} else {
				$stats->{'qiom_obs'} = 100;
			}
		}

		$stats_num_prd_iomsegment_in_this_method = $stats_num_prd_tms_in_this_method + $stats_num_prd_iosegment_in_this_method;
		$stats_num_correct_prd_iomsegment_in_this_method = $stats_num_correct_prd_tms_in_this_method + $stats_num_correct_prd_iosegment_in_this_method;
		if ($stats_num_prd_iomsegment_in_this_method > 0) {
			$stats->{'qiom_prd'} = $stats_num_correct_prd_iomsegment_in_this_method * 100 / $stats_num_prd_iomsegment_in_this_method;
			$stats->{'qiom_prd'} = sprintf("%0.f", $stats->{'qiom_prd'});
		} else {
			if ($stats_num_correct_prd_iomsegment_in_this_method == 0) {
				$stats->{'qiom_prd'} = 0;
			} else {
				$stats->{'qiom_prd'} = 100;
			}
		}

		if ($stats_num_obs_iosegment_in_this_method > 0) {
			$stats->{'qio_obs'} = $stats_num_correct_obs_iosegment_in_this_method * 100 / $stats_num_obs_iosegment_in_this_method;
			$stats->{'qio_obs'} = sprintf("%0.f", $stats->{'qio_obs'});
		} else {
			if ($stats_num_correct_obs_iosegment_in_this_method == 0) {
				$stats->{'qio_obs'} = 0;
			} else {
				$stats->{'qio_obs'} = 100;
			}
		}

		if ($stats_num_prd_iosegment_in_this_method > 0) {
			$stats->{'qio_prd'} = $stats_num_correct_prd_iosegment_in_this_method * 100 / $stats_num_prd_iosegment_in_this_method;
			$stats->{'qio_prd'} = sprintf("%0.f", $stats->{'qio_prd'});
		} else {
			if ($stats_num_correct_prd_iosegment_in_this_method == 0) {
				$stats->{'qio_prd'} = 0;
			} else {
				$stats->{'qio_prd'} = 100;
			}
		}

		if ($stats_num_seqs > 0 ) {
			$stats->{'q3'} = $stats_ratio_correct_iomtopology_residues_in_this_method * 100 / $stats_num_seqs;
			$stats->{'q3'} = sprintf("%0.f", $stats->{'q3'});
		} else {
			$stats->{'q3'} = 0;
		}

		if ($stats_num_seqs > 0 ) {
			$stats->{'q2_iores'} = $stats_ratio_correct_iotopology_residues_in_this_method * 100 / $stats_num_seqs;
			$stats->{'q2_iores'} = sprintf("%0.f", $stats->{'q2_iores'});
		} else {
			$stats->{'q2_iores'} = 0;
		}

		my $io_mcc_TP = $stats_num_correct_prd_i_topology_residues_in_this_method;
		my $io_mcc_TN = $stats_num_correct_prd_o_topology_residues_in_this_method;
		my $io_mcc_FP = $stats_num_i_topology_residues_in_this_method - $stats_num_correct_prd_i_topology_residues_in_this_method;
		my $io_mcc_FN = $stats_num_o_topology_residues_in_this_method - $stats_num_correct_prd_o_topology_residues_in_this_method;
		my $io_mcc_top = ($io_mcc_TP * $io_mcc_TN) - ($io_mcc_FP * $io_mcc_FN);
		my $io_mcc_bottom = sqrt( ($io_mcc_TP + $io_mcc_FP) * ($io_mcc_TP + $io_mcc_FN) * ($io_mcc_TN + $io_mcc_FP) * ($io_mcc_TN + $io_mcc_FN) );
		if ($io_mcc_bottom != 0) {
			$stats->{'io_mcc'} = $io_mcc_top / $io_mcc_bottom;
		}
	}

	return $stats;
}

#===================================================================================================
#	get_aaseq_for_output
#===================================================================================================

sub get_aaseq_for_output {

	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel', 'soluble' );
	my %output_aaseq; # will contain all possible ids as the key and their sequences
	$output->{'aaseq'} = \%output_aaseq;

	#===================================================================================================
	# get the amino acid residue sequence for each sequence in the requested set of protein-chain ids
	#===================================================================================================

	# 3h1j_J:::::___ALLRQAYSALFRRTSTFALTVVLGAVLFERAFDQGADAIFEHLNEGKLWKHIKHKYEASEE:::::

	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {
			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'pdb_seq'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
						my $this_aaseq = trim($bits[1]);
						$output_aaseq{$this_id} = $this_aaseq;
					}
				}
			}
		}
	}
}

#===================================================================================================
#	get_names_for_output
#===================================================================================================

sub get_names_for_output {

	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel', 'soluble' );
	my %output_name; # will contain all possible ids as the key and their names
	$output->{'name'} = \%output_name;

	#===================================================================================================
	# get the names for each sequence in the requested set of protein-chain ids
	#===================================================================================================

	# 3mk7_A:::::The structure of CBB3 cytochrome oxidase:::::Cytochrome c oxidase, cbb3-type, subunit N:::::
	# 3mk7_B:::::The structure of CBB3 cytochrome oxidase:::::Cytochrome c oxidase, cbb3-type, subunit O:::::
	# 3mk7_C:::::The structure of CBB3 cytochrome oxidase:::::Cytochrome c oxidase, cbb3-type, subunit P:::::

	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {
			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'pdbe_name'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
						my $this_name = trim($bits[1]);
						$output_name{$this_id} = $this_name;
					}
				}
			}
		}
	}
}

#===================================================================================================
#	get_release_dates_for_output
#===================================================================================================

sub get_release_dates_for_output {

	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel', 'soluble' );
	my %output_release_dates; # will contain all possible ids as the key and their PDB model release-dates
	$output->{'release_dates'} = \%output_release_dates;

	#===================================================================================================
	# get the PDB model release-dates for each sequence in the requested set of protein-chain ids
	#===================================================================================================

	# 3pjs_K:::::2011-07-06:::::

	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {
			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'release_date'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
						my $this_release_date = trim($bits[1]);
						$output_release_dates{$this_id} = $this_release_date;
					}
				}
			}
		}
	}
}

#===================================================================================================
#	get_methods_resolutions_for_output
#===================================================================================================

sub get_methods_resolutions_for_output {

	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel', 'soluble' );
	my %output_pdb_methods; # will contain all possible ids as the key and their PDB model info
	my %output_pdb_resolutions; # will contain all possible ids as the key and their PDB model info
	$output->{'pdb_methods'} = \%output_pdb_methods;
	$output->{'pdb_resolutions'} = \%output_pdb_resolutions;

	#===================================================================================================
	# get the PDB model release-dates for each sequence in the requested set of protein-chain ids
	#===================================================================================================

	# 3pjs_K:::::2011-07-06:::::

	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {
			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'pdb_method_resolution'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					# 1mhs_A:::::METHOD,,,,,ELECTRON CRYSTALLOGRAPHY;;;;;RESOLUTION,,,,,8.0;;;;;:::::
					# 1n7l_A:::::METHOD,,,,,SOLUTION NMR;;;;;:::::
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists ($output->{'id'}->{$this_protein_category}->{$this_id}) ) {
						my @bits2 = split(/\;\;\;\;\;/, $bits[1]);
						foreach my $this_method_resolution (@bits2) {
							my @bits2 = split(/\,\,\,\,\,/, $this_method_resolution);
							my $this_name = $bits2[0];
							my $this_value = $bits2[1];
							if ($this_name eq 'METHOD') {
								$output_pdb_methods{$this_id} = $this_value;
							} elsif ($this_name = 'RESOLUTION') {
								$output_pdb_resolutions{$this_id} = $this_value;
							}
						}
					}
				}
			}
		}
	}
}

#===================================================================================================
#	get_mpstruc_category_for_output
#===================================================================================================

sub get_mpstruc_category_for_output {

	my %output_mpstruc_category;
	$output->{'mpstruc_category'} = \%output_mpstruc_category;

	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel', 'soluble' );
	my @protein_category_names = ( 'Alpha-Helical Membrane Protein', 'Beta-Barrel Membrane Protein', 'Soluble Protein' );

	my $i = 0;
	for my $this_protein_category ( @protein_categories ) {
		for my $this_id ( %{$output->{'id'}->{$this_protein_category}} ) {
			$output_mpstruc_category{$this_id} = '';
			if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
				$output_mpstruc_category{$this_id} = $protein_category_names[$i];
			}
		}
		$i++;
	}

	#===================================================================================================
	# get the category of each protein as assigned by the White lab's Membrane Proteins of Known Structure
	#===================================================================================================

	# 1ehk_A:::::Electron Transport Chain Complexes: Complex IV (Cytochrome C Oxidase):::::

	if ($global->{'kingdom'}->{'alpha_polytopic_bitopic'} ne '') {
		my $file_having_list_of_ids = $data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'mpstruc_category'};
		open LISTFILE, $file_having_list_of_ids or die $!;
		my @input_lines = <LISTFILE>;
		foreach my $input_line (@input_lines) {
			chomp $input_line;
			if ($input_line ne '') {
				my @bits = split(/:::::/, $input_line);
				my $this_id = trim($bits[0]);
				if ( exists $output->{'id'}->{'alpha_polytopic_bitopic'}->{$this_id} ) {
					my $this_mpstruc_category = trim($bits[1]);
					$output_mpstruc_category{$this_id} = $this_mpstruc_category;
				}
			}
		}
	}
}

#===================================================================================================
#	get_opm_tm_segments_for_output
#===================================================================================================

sub get_opm_tm_segments_for_output {

	my $run_type = shift;
	my %output_opm_tm_segments; # will contain all possible ids as the key and their OPM segments
	my %output_opm_tm_segments_topology; # will contain all possible ids as the key and their OPM segments
	$output->{'opm_tm_segments'} = \%output_opm_tm_segments;
	$output->{'opm_tm_segments_topology'} = \%output_opm_tm_segments_topology;

	for my $this_id ( %{$output->{'id'}->{'all'}} ) {
		$output_opm_tm_segments{$this_id} = '';
		$output_opm_tm_segments_topology{$this_id} = '';
	}

	#===================================================================================================
	# get the OPM transmembrane helix (TMH) and transmembrane beta-barrel (TMBB) segments for each sequence in the requested set of protein-chain ids
	#===================================================================================================

	# 1fjk_A:::::____________________________HHHHHHHHHHHHHHHHHHHHHHH_:::::
	# 1grm_A:::::__BBBBBBBBBBBBB_:::::

	my @protein_categories;
	if ($run_type eq 'benchmark') {
		@protein_categories = ( 'alpha_polytopic_bitopic' );
	} elsif ($run_type eq 'view') {
		@protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel' )
	}
	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {

			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'opm_tm_segments'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
						my $this_aaseq = trim($bits[1]);
						$output_opm_tm_segments{$this_id} = $this_aaseq;
					}
				}
			}

			$file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'opm_tm_segments_topology'};
			open LISTFILE2, $file_having_list_of_ids or die $!;
			my @input_lines2 = <LISTFILE2>;
			foreach my $input_line (@input_lines2) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
						my $this_aaseq = trim($bits[1]);
						$output_opm_tm_segments_topology{$this_id} = $this_aaseq;
					}
				}
			}
		}
	}
}

#===================================================================================================
#	get_pdbtm_segments_for_output
#===================================================================================================

sub get_pdbtm_segments_for_output {

	my $run_type = shift;
	my %output_pdbtm_segments; # will contain all possible ids as the key and their PDBTM segments
	my %output_pdbtm_segments_topology; # will contain all possible ids as the key and their PDBTM segments
	$output->{'pdbtm_segments'} = \%output_pdbtm_segments;
	$output->{'pdbtm_segments_topology'} = \%output_pdbtm_segments_topology;
	my %output_aaseq = %{$output->{'aaseq'}};

	for my $this_id ( %{$output->{'id'}->{'all'}} ) {
		$output_pdbtm_segments{$this_id} = '';
		$output_pdbtm_segments_topology{$this_id} = '';
	}

	#===================================================================================================
	# get the PDBTM segments for each sequence in the requested set of protein-chain ids
	#===================================================================================================

	# one line (split here for viewing) :
	# 1h6i_A:::::ALPHA__side_1,,,,,9-13,69-71,87-93,157-167,232-233,;;;;;ALPHA__side_2,,,,,33-49,115-138,185-187,203-209,;;;;;
	# ALPHA__alpha_helix,,,,,14-32,50-68,94-114,139-156,168-184,210-231,;;;;;ALPHA__membrane_loop,,,,,72-86,188-202,;;;;;ALPHA__unknown,,,,,1-8,234-269,;;;;;:::::

	# ALPHA__side_1,,,,,47-47,93-97,145-163,168-177,194-194,230-238,286-290,
	# ALPHA__side_2,,,,,-1-32,61-79,115-129,206-206,207-211,252-265,305-334,
	# ALPHA__alpha_helix,,,,,33-46,48-60,80-92,98-114,130-144,195-205,212-229,239-251,266-285,291-304,
	# ALPHA__membrane_loop,,,,,178-193,
	# ALPHA__unknown,,,,,-2--2,164-167,335-337,

	#my %type_code = (ALPHA__side_1 => "1", ALPHA__side_2 => "2", ALPHA__alpha_helix => "H", ALPHA__membrane_loop => "L", ALPHA__coil => "C", ALPHA__membrane_inside => "I", ALPHA__beta_strand => "B", ALPHA__unknown => "U" );
	my %type_code = (side_1 => "1", side_2 => "2", alpha_helix => "H", membrane_loop => "L", coil => "C", membrane_inside => "I", beta_strand => "B", unknown => "U" );

	my @protein_categories;
	if ($run_type eq 'benchmark') {
		@protein_categories = ( 'alpha_polytopic_bitopic' );
	} elsif ($run_type eq 'view') {
		@protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel' )
	}
	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {

			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'pdbtm_segments'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {

						my $this_pdbtm_segments = trim($bits[1]);
						my $length_aaseq = length($output_aaseq{$this_id});
						#my $this_pdbtm_segments_by_residue = ( '_' x $length_aaseq );
						my @this_pdbtm_segments_by_residue_array;
						for ( my $j = 0; $j < $length_aaseq; $j++ ) {
							$this_pdbtm_segments_by_residue_array[$j] = '_';
						}

						if ($this_pdbtm_segments ne '-') {
							my @bits2 = split(/;;;;;/, $this_pdbtm_segments);
							foreach my $bit2 (@bits2) {
								if ($bit2 ne '') {
									my @bits3 = split(/,,,,,/, $bit2);
									my $this_type = $bits3[0];
									my $this_values = $bits3[1];
									my @bits3a = split(/__/, $this_type);
									$this_type = $bits3a[1];
									my $this_char = '_';
									if ( exists($type_code{$this_type}) ) {
										$this_char = $type_code{$this_type};
									}
									my @bits4 = split(/,/, $this_values);
									foreach my $bit4 (@bits4) {
										$bit4 =~ m/(-?\d+)-(-?\d+)/;
										my $from = trim($1);
										my $to = trim($2);
										for ( my $i = ($from - 1); $i < $to; $i++ ) {
											if ($i >= 0) {
												$this_pdbtm_segments_by_residue_array[$i] = $this_char;
											}
										}
									}
								}
							}
						}

						my $this_pdbtm_segments_by_residue = '';
						for ( my $j = 0; $j < $length_aaseq; $j++ ) {
							$this_pdbtm_segments_by_residue .= $this_pdbtm_segments_by_residue_array[$j];
						}
						$output_pdbtm_segments{$this_id} = $this_pdbtm_segments_by_residue;
					}
				}
			}

			$file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'pdbtm_segments_topology'};
			open LISTFILE2, $file_having_list_of_ids or die $!;
			my @input_lines2 = <LISTFILE2>;
			foreach my $input_line (@input_lines2) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
						my $this_aaseq = trim($bits[1]);
						$output_pdbtm_segments_topology{$this_id} = $this_aaseq;
					}
				}
			}
		}
	}
}

#===================================================================================================
#	get_opm_adj_tm_segments_for_output
#===================================================================================================

sub get_opm_adj_tm_segments_for_output {

	my $run_type = shift;
	my %output_opm_adj_tm_segments; # will contain all possible ids as the key and their OPM segments
	my %output_opm_adj_tm_segments_topology; # will contain all possible ids as the key and their OPM segments
	$output->{'opm_adj_tm_segments'} = \%output_opm_adj_tm_segments;
	$output->{'opm_adj_tm_segments_topology'} = \%output_opm_adj_tm_segments_topology;

	for my $this_id ( %{$output->{'id'}->{'all'}} ) {
		$output_opm_adj_tm_segments{$this_id} = '';
		$output_opm_adj_tm_segments_topology{$this_id} = '';
	}

	#===================================================================================================
	# get the OPM transmembrane helix (TMH) and transmembrane beta-barrel (TMBB) segments for each sequence in the requested set of protein-chain ids
	#===================================================================================================

	# 1fjk_A:::::____________________________HHHHHHHHHHHHHHHHHHHHHHH_:::::
	# 1grm_A:::::__BBBBBBBBBBBBB_:::::

	my @protein_categories;
	if ($run_type eq 'benchmark') {
		@protein_categories = ( 'alpha_polytopic_bitopic' );
	} elsif ($run_type eq 'view') {
		@protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel' )
	}
	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {

			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'opm_adj_tm_segments'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
						my $this_aaseq = trim($bits[1]);
						$output_opm_adj_tm_segments{$this_id} = $this_aaseq;
					}
				}
			}

			$file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'opm_adj_tm_segments_topology'};
			open LISTFILE2, $file_having_list_of_ids or die $!;
			my @input_lines2 = <LISTFILE2>;
			foreach my $input_line (@input_lines2) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
						my $this_aaseq = trim($bits[1]);
						$output_opm_adj_tm_segments_topology{$this_id} = $this_aaseq;
					}
				}
			}
		}
	}
}

#===================================================================================================
#	get_PDB_DSSP_for_output
#===================================================================================================

sub get_PDB_DSSP_for_output {

	my %output_dssp; # will contain all possible ids as the key and their sequences
	$output->{'dssp'} = \%output_dssp;

	for my $this_id ( %{$output->{'id'}->{'all'}} ) {
		$output_dssp{$this_id} = '';
	}

	#===================================================================================================
	# get the DSSP secondary structure from PDB website for each sequence in the requested set of protein-chain ids
	#===================================================================================================

	# 1jb0_L:::::_____EEGGG_TTBB_EE_HHHH_HHHHHHHHTSTTT_TT__HHHHHHHHHHHHHHHHTHHHHHHSTTTTSTTHHHHHHHHHHHHHHHHHHHHHHHHHHHTSS___SS_GGGSHHHHHHHHHHHHHHHHHHHHHHHHHHHTHHHHHHHHHTTT_:::::
	# 1jb0_X:::::__________HHHHHHHHHHHHHHHHHHHHHTS__:::::

	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel', 'soluble' );

	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {

			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'dssp'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
						my $this_aaseq = trim($bits[1]);
						$output_dssp{$this_id} = $this_aaseq;
					}
				}
			}
		}
	}
}

#===================================================================================================
#	get_PDB_STRIDE_for_output
#===================================================================================================

sub get_PDB_STRIDE_for_output {

	my %output_stride; # will contain all possible ids as the key and their sequences
	$output->{'stride'} = \%output_stride;

	for my $this_id ( %{$output->{'id'}->{'all'}} ) {
		$output_stride{$this_id} = '';
	}

	#===================================================================================================
	# get the STRIDE secondary structure that the STRIDE program read from the OPM PDB files
	#===================================================================================================

	# 1jb0_L:::::_____EEGGG_TTBB_EE_HHHH_HHHHHHHHTSTTT_TT__HHHHHHHHHHHHHHHHTHHHHHHSTTTTSTTHHHHHHHHHHHHHHHHHHHHHHHHHHHTSS___SS_GGGSHHHHHHHHHHHHHHHHHHHHHHHHHHHTHHHHHHHHHTTT_:::::
	# 1jb0_X:::::__________HHHHHHHHHHHHHHHHHHHHHTS__:::::

	#my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel', 'soluble' );
	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel' );
	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {

			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'stride'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
						my $this_aaseq = trim($bits[1]);
						$output_stride{$this_id} = $this_aaseq;
					}
				}
			}
		}
	}
}

#===================================================================================================
#	get_PDB_file_seq_for_output
#===================================================================================================

sub get_PDB_file_seq_for_output {

	my %output_aaseq; # will contain all possible ids as the key and their sequences
	$output->{'opm_seq'} = \%output_aaseq;

	for my $this_id ( %{$output->{'id'}->{'all'}} ) {
		$output_aaseq{$this_id} = '';
	}

	#===================================================================================================
	# get the amino acid sequence from OPM's PDB file for each sequence in the requested set of protein-chain ids
	#===================================================================================================

	# 1zoy_D:::::_________________________________ASSKAASLHWTGERVVSVLLLGLLPAAYLNPCSAMDYSLAAALTLHGHWGIGQVVTDYVRGDALQKAAKAGLLALSAFTFAGLCYFNYHDVGICKAVAMLWKL:::::

	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel', 'soluble' );
	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {

			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'opm_seq'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
						my $this_aaseq = trim($bits[1]);
						$output_aaseq{$this_id} = $this_aaseq;
					}
				}
			}
		}
	}
}

#===================================================================================================
#	get_PDB_OPM_STRIDE_file_seq_for_output
#===================================================================================================

sub get_PDB_OPM_STRIDE_file_seq_for_output {

	my %output_aaseq; # will contain all possible ids as the key and their sequences
	$output->{'stride_seq'} = \%output_aaseq;

	for my $this_id ( %{$output->{'id'}->{'all'}} ) {
		$output_aaseq{$this_id} = '';
	}

	#===================================================================================================
	# get the amino acid sequence from OPM's PDB file for each sequence as was read by STRIDE
	#===================================================================================================

	# 1zoy_D:::::_________________________________ASSKAASLHWTGERVVSVLLLGLLPAAYLNPCSAMDYSLAAALTLHGHWGIGQVVTDYVRGDALQKAAKAGLLALSAFTFAGLCYFNYHDVGICKAVAMLWKL:::::

	#my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel', 'soluble' );
	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel' );
	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {

			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'stride_seq'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
						my $this_aaseq = trim($bits[1]);
						$output_aaseq{$this_id} = $this_aaseq;
					}
				}
			}
		}
	}
}

#===================================================================================================
#	get_OPM_tmh_tmbb_for_output
#===================================================================================================

sub get_OPM_tmh_tmbb_for_output {

	my %output_opm_tmh_tmbb_segments; # will contain all possible ids as the key and their sequences
	$output->{'opm_tmh_tmbb_segments'} = \%output_opm_tmh_tmbb_segments;

	for my $this_id ( %{$output->{'id'}->{'all'}} ) {
		$output_opm_tmh_tmbb_segments{$this_id} = '';
	}

	#===================================================================================================
	# get the OPM membrane segments that have been set to HHHHHH or BBBBBB by comparing to DSSP from PDB website
	#===================================================================================================

	# 1a11_A:::::iiiiHHHHHHHHHHHHHHHHHHooo:::::

	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel' );
	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {

			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'opm_tmh_tmbb_segments'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
						my $this_aaseq = trim($bits[1]);
						$output_opm_tmh_tmbb_segments{$this_id} = $this_aaseq;
					}
				}
			}
		}
	}
}

#===================================================================================================
#	get_OPM_adj_tmh_tmbb_for_output
#===================================================================================================

sub get_OPM_adj_tmh_tmbb_for_output {

	my %output_opm_adj_tmh_tmbb_segments; # will contain all possible ids as the key and their sequences
	$output->{'opm_adj_tmh_tmbb_segments'} = \%output_opm_adj_tmh_tmbb_segments;

	for my $this_id ( %{$output->{'id'}->{'all'}} ) {
		$output_opm_adj_tmh_tmbb_segments{$this_id} = '';
	}

	#===================================================================================================
	# get the OPM membrane segments that have been set to HHHHHH or BBBBBB by comparing to DSSP from PDB website
	#===================================================================================================

	# 1a11_A:::::iiiiHHHHHHHHHHHHHHHHHHooo:::::

	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel' );
	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {

			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'opm_adj_tmh_tmbb_segments'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
						my $this_aaseq = trim($bits[1]);
						$output_opm_adj_tmh_tmbb_segments{$this_id} = $this_aaseq;
					}
				}
			}
		}
	}
}

#===================================================================================================
#	get_OPM_tm_subunits_for_output
#===================================================================================================

sub get_OPM_tm_subunits_for_output {

	my %output_opm_segments; # will contain all possible ids as the key and their sequences
	$output->{'opm_tm_subunits'} = \%output_opm_segments;
	my %output_aaseq = %{$output->{'aaseq'}};

	for my $this_id ( %{$output->{'id'}->{'all'}} ) {
		$output_opm_segments{$this_id} = '';
	}

	#===================================================================================================
	# get the OPM transmembrane segments (includes helices and beta-sheets) for each sequence in the requested set of protein-chain ids
	#===================================================================================================

	# 3eml_A:::::OPM-TRANSMEMBRANE-SUBUNIT,,,,,6-29,42-65,79-100,120-141,176-198,234-255,268-290,;;;;;:::::
	# 3f5w_C:::::OPM-TRANSMEMBRANE-SUBUNIT,,,,,31-50,86-110,;;;;;:::::

	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel' );
	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {

			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'opm_tm_subunits'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {

						my $this_opm_segments = trim($bits[1]);
						my $length_aaseq = length($output_aaseq{$this_id});
						#my $this_opm_segments_by_residue = ( '_' x $length_aaseq );
						my @this_opm_segments_by_residue_array;
						for ( my $j = 0; $j < $length_aaseq; $j++ ) {
							$this_opm_segments_by_residue_array[$j] = '_';
						}

						if ($this_opm_segments ne '-') {
							my @bits2 = split(/;;;;;/, $this_opm_segments);
							foreach my $bit2 (@bits2) {
								if ($bit2 ne '') {
									my @bits3 = split(/,,,,,/, $bit2);
									my $this_type = $bits3[0];
									my $this_values = $bits3[1];
									if ($this_type eq 'OPM-TRANSMEMBRANE-SUBUNIT') {
										my $this_char = 'm';
										my @bits4 = split(/,/, $this_values);
										foreach my $bit4 (@bits4) {
											$bit4 =~ m/(-?\d+)-(-?\d+)/;
											my $from = trim($1);
											my $to = trim($2);
											for ( my $i = ($from - 1); $i < $to; $i++ ) {
												if ($i >= 0) {
													$this_opm_segments_by_residue_array[$i] = $this_char;
												}
											}
										}
									}
								}
							}
						}

						my $this_opm_segments_by_residue = '';
						for ( my $j = 0; $j < $length_aaseq; $j++ ) {
							$this_opm_segments_by_residue .= $this_opm_segments_by_residue_array[$j];
						}
						$output_opm_segments{$this_id} = $this_opm_segments_by_residue;
					}
				}
			}
		}
	}
}

#===================================================================================================
#	get_OPM_topology_for_output
#===================================================================================================

sub get_OPM_topology_for_output {

	my %output_opm_segments; # will contain all possible ids as the key and their OPM topology by residue
	$output->{'opm_posn'} = \%output_opm_segments;
	my %output_aaseq = %{$output->{'aaseq'}};

	for my $this_id ( %{$output->{'id'}->{'all'}} ) {
		$output_opm_segments{$this_id} = '';
	}

	#===================================================================================================
	# get the OPM-assigned topology for each sequence in the requested set of protein-chain ids
	#===================================================================================================

	# 1a11_A:::::MEMBRANE-SEGMENTS,,,,,5-22,;;;;;INSIDE-SEGMENTS,,,,,1-4,;;;;;OUTSIDE-SEGMENTS,,,,,23-25,;;;;;:::::
	# 1a91_A:::::MEMBRANE-SEGMENTS,,,,,8-34,49-49,51-75,;;;;;INSIDE-SEGMENTS,,,,,1-7,76-79,;;;;;OUTSIDE-SEGMENTS,,,,,35-48,50-50,;;;;;:::::

	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel' );
	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {

			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'opm_posn'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {

						my $this_opm_segments = trim($bits[1]);
						my $length_aaseq = length($output_aaseq{$this_id});
						#my $this_opm_segments_by_residue = ( '_' x $length_aaseq );
						my @this_opm_segments_by_residue_array;
						for ( my $j = 0; $j < $length_aaseq; $j++ ) {
							$this_opm_segments_by_residue_array[$j] = '_';
						}

						if ($this_opm_segments ne '-') {
							my @bits2 = split(/;;;;;/, $this_opm_segments);
							foreach my $bit2 (@bits2) {
								if ($bit2 ne '') {
									my @bits3 = split(/,,,,,/, $bit2);
									my $this_type = $bits3[0];
									my $this_values = $bits3[1];
									my $this_char = '?';
									if ($this_type eq 'MEMBRANE-SEGMENTS') {
										$this_char = 'm';
									} elsif ($this_type eq 'INSIDE-SEGMENTS') {
										$this_char = 'i';
									} elsif ($this_type eq 'OUTSIDE-SEGMENTS') {
										$this_char = 'o';
									}
									my @bits4 = split(/,/, $this_values);
									foreach my $bit4 (@bits4) {
										$bit4 =~ m/(-?\d+)-(-?\d+)/;
										my $from = trim($1);
										my $to = trim($2);
										for ( my $i = ($from - 1); $i < $to; $i++ ) {
											if ($i >= 0) {
												$this_opm_segments_by_residue_array[$i] = $this_char;
											}
										}
									}
								}
							}
						}

						my $this_opm_segments_by_residue = '';
						for ( my $j = 0; $j < $length_aaseq; $j++ ) {
							$this_opm_segments_by_residue .= $this_opm_segments_by_residue_array[$j];
						}
						$output_opm_segments{$this_id} = $this_opm_segments_by_residue;
					}
				}
			}
		}
	}
}

#===================================================================================================
#	get_OPM_membrane_near_for_output
#===================================================================================================

sub get_OPM_membrane_near_for_output {

	my %output_aaseq; # will contain all possible ids as the key and their sequences
	$output->{'opm_membrane_near'} = \%output_aaseq;

	for my $this_id ( %{$output->{'id'}->{'all'}} ) {
		$output_aaseq{$this_id} = '';
	}

	#===================================================================================================
	# get the amino acid sequence from OPM's PDB file for each sequence in the requested set of protein-chain ids
	#===================================================================================================

	# 1a11_A:::::____IIIIIIIIIIOOOOOOOO___:::::

	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel' );
	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {

			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'opm_membrane_near'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
						my $this_aaseq = trim($bits[1]);
						$output_aaseq{$this_id} = $this_aaseq;
					}
				}
			}
		}
	}
}

#===================================================================================================
#	get_OPM_tm_segments_topology_for_output
#===================================================================================================

sub get_OPM_tm_segments_topology_for_output {

	my %output_aaseq; # will contain all possible ids as the key and their sequences
	$output->{'opm_tm_segments_topology'} = \%output_aaseq;

	for my $this_id ( %{$output->{'id'}->{'all'}} ) {
		$output_aaseq{$this_id} = '';
	}

	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel' );
	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {

			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'opm_tm_segments_topology'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
						my $this_aaseq = trim($bits[1]);
						$output_aaseq{$this_id} = $this_aaseq;
					}
				}
			}
		}
	}
}

#===================================================================================================
#	get_OPM_adj_tm_segments_topology_for_output
#===================================================================================================

sub get_OPM_adj_tm_segments_topology_for_output {

	my %output_aaseq; # will contain all possible ids as the key and their sequences
	$output->{'opm_adj_tm_segments_topology'} = \%output_aaseq;

	for my $this_id ( %{$output->{'id'}->{'all'}} ) {
		$output_aaseq{$this_id} = '';
	};

	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel' );
	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} ne '') {

			my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'opm_adj_tm_segments_topology'};
			open LISTFILE, $file_having_list_of_ids or die $!;
			my @input_lines = <LISTFILE>;
			foreach my $input_line (@input_lines) {
				chomp $input_line;
				if ($input_line ne '') {
					my @bits = split(/:::::/, $input_line);
					my $this_id = trim($bits[0]);
					if ( exists $output->{'id'}->{$this_protein_category}->{$this_id} ) {
						my $this_aaseq = trim($bits[1]);
						$output_aaseq{$this_id} = $this_aaseq;
					}
				}
			}
		}
	}
}

#===================================================================================================
#	read_and_validate_input_parameters
#===================================================================================================

sub read_and_validate_input_parameters {

	#===================================================================================================
	# read the input parameters from the user
	#===================================================================================================

	my $params = $q->Vars;

	my $benchmark_standard = '';
	$global->{'minimum_overlap'} = $data->{'minimum_overlap'};
	$global->{'minimum_helix_length'} = $data->{'minimum_helix_length'};
	$global->{'max_dist_approx_ends'} = $data->{'max_dist_approx_ends'};
	$global->{'display_width'} = $data->{'display_width'};
	$global->{'sort_by'} = $data->{'sort_by'};
	$global->{'sort_by_name'} = $data->{'sort_by_name'};
	$global->{'graph_by'} = '';
	$global->{'graph_by_2'} = '';
	$global->{'graph_by_3'} = '';
	$global->{'graph_title'} = '';
	$global->{'tmh_numTMH'} = '';
	$global->{'helix_type'} = '';
	$global->{'tmh_family'} = '';
	$global->{'num_soluble'} = 'all';
	$global->{'num_output_display_methods'} = 0;
	$global->{'num_output_prediction_methods'} = 0;
	$global->{'num_output_sequences'} = 0;
	my $manual_list_seqs = '';
	my $run = '';

	my %global_similarity_percent;
	$global->{'similarity_percent'} = \%global_similarity_percent;
	$global->{'similarity_percent'}->{'alpha_polytopic_bitopic'} = '';
	$global->{'similarity_percent'}->{'beta_barrel'} = '';
	$global->{'similarity_percent'}->{'soluble'} = '';

	my %global_water_needle;
	$global->{'water_needle'} = \%global_water_needle;
	$global->{'water_needle'}->{'alpha_polytopic_bitopic'} = '';
	$global->{'water_needle'}->{'beta_barrel'} = '';
	$global->{'water_needle'}->{'soluble'} = '';

	my %global_similarity_identity;
	$global->{'similarity_identity'} = \%global_similarity_identity;
	$global->{'similarity_identity'}->{'alpha_polytopic_bitopic'} = '';
	$global->{'similarity_identity'}->{'beta_barrel'} = '';
	$global->{'similarity_identity'}->{'soluble'} = '';

	my %global_kingdom;
	$global->{'kingdom'} = \%global_kingdom;
	$global->{'kingdom'}->{'alpha_polytopic_bitopic'} = '';
	$global->{'kingdom'}->{'beta_barrel'} = '';
	$global->{'kingdom'}->{'soluble'} = '';

	my %global_nonhomologous_list;
	$global->{'nonhomologous_list'} = \%global_nonhomologous_list;
	$global->{'nonhomologous_list'}->{'alpha_polytopic_bitopic'} = '';
	$global->{'nonhomologous_list'}->{'beta_barrel'} = '';
	$global->{'nonhomologous_list'}->{'soluble'} = '';

	$global->{'restrict_experimental_method'} = '';
	$global->{'restrict_resolution'} = '';

	if (defined($params->{'benchmark_standard'})) {
		$global->{'benchmark_standard'} = trim( $params->{'benchmark_standard'} );
	}
	if (defined($params->{'tmh_similarity_percent'})) {
		$global->{'similarity_percent'}->{'alpha_polytopic_bitopic'} = trim( $params->{'tmh_similarity_percent'} );
	}
	if (defined($params->{'tmh_water_needle'})) {
		$global->{'water_needle'}->{'alpha_polytopic_bitopic'} = trim( $params->{'tmh_water_needle'} );
	}
	if (defined($params->{'tmh_similarity_identity'})) {
		$global->{'similarity_identity'}->{'alpha_polytopic_bitopic'} = trim( $params->{'tmh_similarity_identity'} );
	}
	if (defined($params->{'tmh_kingdom'})) {
		$global->{'kingdom'}->{'alpha_polytopic_bitopic'} = trim( $params->{'tmh_kingdom'} );
		if ( lc($global->{'kingdom'}->{'alpha_polytopic_bitopic'}) eq 'none' ) {
			$global->{'kingdom'}->{'alpha_polytopic_bitopic'} = '';
		}
	}
	if (defined($params->{'tmh_numTMH'})) {
		$global->{'tmh_numTMH'} = trim( $params->{'tmh_numTMH'} );
	}
	if (defined($params->{'helix_type'})) {
		$global->{'helix_type'} = trim( $params->{'helix_type'} );
	}
	if (defined($params->{'tmh_family'})) {
		$global->{'tmh_family'} = trim( $params->{'tmh_family'} );
	}
	if (defined($params->{'restrict_year'})) {
		$global->{'restrict_year'} = trim( $params->{'restrict_year'} );
		if ($global->{'restrict_year'} eq 'none') {
			$global->{'restrict_year'} = '';
		}
	}
	if (defined($params->{'restrict_resolution'})) {
		$global->{'restrict_resolution'} = trim( $params->{'restrict_resolution'} );
		if ( (is_valid_positive_number_from_zero($global->{'restrict_resolution'})) == 0 ) {
			$global->{'got_error'} = 1;
			$global->{'err_msg'} .= "The requested resolution limit : " . $global->{'restrict_resolution'} . " is not a positive number.<br>";
		}
	}
	if ( (defined($params->{'xray_diffraction'})) || (defined($params->{'solution_nmr'})) || (defined($params->{'electron_crystallography'})) || (defined($params->{'electron_microscopy'})) || (defined($params->{'fiber_diffraction'})) || (defined($params->{'solid_state_nmr'})) ) {
		$global->{'restrict_experimental_method'} = ';;;;;';
		if (defined($params->{'xray_diffraction'})) {
			$global->{'restrict_experimental_method'} .= 'X-RAY DIFFRACTION;;;;;';
		}
		if (defined($params->{'solution_nmr'})) {
			$global->{'restrict_experimental_method'} .= 'SOLUTION NMR;;;;;';
		}
		if (defined($params->{'electron_crystallography'})) {
			$global->{'restrict_experimental_method'} .= 'ELECTRON CRYSTALLOGRAPHY;;;;;';
		}
		if (defined($params->{'electron_microscopy'})) {
			$global->{'restrict_experimental_method'} .= 'ELECTRON MICROSCOPY;;;;;';
		}
		if (defined($params->{'fiber_diffraction'})) {
			$global->{'restrict_experimental_method'} .= 'FIBER DIFFRACTION;;;;;';
		}
		if (defined($params->{'solid_state_nmr'})) {
			$global->{'restrict_experimental_method'} .= 'SOLID-STATE NMR;;;;;';
		}
	}
	if (defined($params->{'bb_similarity_percent'})) {
		$global->{'similarity_percent'}->{'beta_barrel'} = trim( $params->{'bb_similarity_percent'} );
	}
	if (defined($params->{'bb_water_needle'})) {
		$global->{'water_needle'}->{'beta_barrel'} = trim( $params->{'bb_water_needle'} );
	}
	if (defined($params->{'bb_similarity_identity'})) {
		$global->{'similarity_identity'}->{'beta_barrel'} = trim( $params->{'bb_similarity_identity'} );
	}
	if (defined($params->{'bb_kingdom'})) {
		$global->{'kingdom'}->{'beta_barrel'} = trim( $params->{'bb_kingdom'} );
		if ( lc($global->{'kingdom'}->{'beta_barrel'}) eq 'none' ) {
			$global->{'kingdom'}->{'beta_barrel'} = '';
		} 
	}
	if (defined($params->{'soluble_kingdom'})) {
		$global->{'kingdom'}->{'soluble'} = lc( trim( $params->{'soluble_kingdom'} ) );
		if ( $global->{'kingdom'}->{'soluble'} eq 'none' ) {
			$global->{'kingdom'}->{'soluble'} = '';
		}
	}
	if (defined($params->{'manual_list_seqs'})) {
		$global->{'manual_list_seqs'} = trim( $params->{'manual_list_seqs'} );
	}
	if (defined($params->{'run'})) {
		$run = trim( $params->{'run'} );
	}

	if ($global->{'got_error'} == 0) {
		if (($run eq 'fasta') or ($run eq 'fastaG') or ($run eq 'list') or ($run eq 'names') or ($run eq 'count')) {
			$global->{'display_text_not_html'} = 1;
			$global->{'run'} = $run;
		} elsif (($run eq 'dload') or ($run eq 'dloadG')) {
			$global->{'display_text_not_html'} = 2;
			$global->{'run'} = $run;
		} elsif (($run eq 'benchmark') or ($run eq 'display') or ($run eq 'display_date')) {
			$global->{'display_text_not_html'} = 0;
			$global->{'run'} = $run;
		}
	}

	my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel', 'soluble' );
	for my $this_protein_category ( @protein_categories ) {
		if ($global->{'kingdom'}->{$this_protein_category} eq '') {
			$global->{$this_protein_category}->{'soluble'} = '';
			$global->{$this_protein_category}->{'soluble'} = '';
			$global->{$this_protein_category}->{'soluble'} = '';
			$global->{$this_protein_category}->{'soluble'} = '';
		}
	}
	if ($global->{'kingdom'}->{'alpha_polytopic_bitopic'} eq '') {
		$global->{'tmh_numTMH'} = '';
		$global->{'helix_type'} = '';
		$global->{'tmh_family'} = '';
		$global->{'restrict_year'} = '';
	}

	if ($global->{'got_error'} == 0) {
		if ($global->{'manual_list_seqs'} ne '') {

			my %entered_id_is_ok;
			my @bits = split( /\n|\r/, $global->{'manual_list_seqs'} );
			foreach my $this_input_id (@bits) {
				if ($this_input_id ne '') {
					$this_input_id = trim($this_input_id);
					my $this_pdbid = substr( $this_input_id, 0, 4 );
					$this_pdbid = lc($this_pdbid);
					my $this_chain = substr( $this_input_id, 4, 2 );
					$this_input_id = $this_pdbid . $this_chain;
					$entered_id_is_ok{$this_input_id} = 0;
				}
			}

			my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel', 'soluble' );
			for my $this_protein_category ( @protein_categories ) {

				my $user_has_sequences_from_this_protein_category = 0;
				my $file_having_list_of_ids = $data->{'file_name'}->{$this_protein_category}->{'id'};
				open LISTFILE, $file_having_list_of_ids or die $!;
				my @input_lines = <LISTFILE>;
				foreach my $input_line (@input_lines) {
					chomp $input_line;
					if ($input_line ne '') {
						my @bits = split(/:::::/, $input_line);
						my $this_reference_id = trim($bits[0]);
						foreach my $entered_id ( keys %entered_id_is_ok ) {
							if ($this_reference_id eq $entered_id) {
								$entered_id_is_ok{$entered_id} = 1;
								$output->{'id'}->{$this_protein_category}->{$entered_id} = '';
								$output->{'id'}->{'all'}->{$entered_id} = '';
								$user_has_sequences_from_this_protein_category = 1;
							}
						}
					}
				}

				# This kingdom field is used as a protein_category flag for validations and calculations,
				# because processing is different for the different protein categories.
				# Benchmarking for alpha-helical membrane proteins requires reading benchmark files,
				# whereas for beta-barrel membrane proteins and soluble proteins it is assumed that there are no membrane helices
				# and so it is not necessary to read benchmark files to get the topography and topology.
				# When the user enters the benchmark list and we find out what protein categories the entered sequences are,
				# then flag those sequences so that the processing of those protein category types will be done.
				if ($user_has_sequences_from_this_protein_category == 1) {
					$global->{'kingdom'}->{$this_protein_category} = 'all';
				}
			}

			my $err_msg = '';
			my $err_list = '';
			foreach my $entered_id ( keys %entered_id_is_ok ) {
				if ($entered_id_is_ok{$entered_id} == 0) {
					$err_msg = 'Entered sequence-id(s) is(are) unknown sequence-id : ';
					$err_list .= $entered_id . ',';
				}
			}
			if ($err_list ne '') {
				chop $err_list;
				if ( $err_list =~ /,/ ) {
					$err_msg = 'Entered sequence-ids are not valid sequence-ids to use : ';
				} else {
					$err_msg = 'Entered sequence-id is not a valid sequence-id to use : ';
				}
				$err_msg = $err_msg . $err_list . '.';
				$global->{'info_msg'} .= $err_msg;
			}
		}
	}

	if (($global->{'got_error'} == 0) and ($global->{'manual_list_seqs'} eq '')) {
		if (($global->{'kingdom'}->{'soluble'} ne 'none') and ($global->{'kingdom'}->{'soluble'} ne '')) {
			if (defined($params->{'num_soluble'})) {
				$global->{'num_soluble'} = lc( trim( $params->{'num_soluble'} ) );
				if ( $global->{'num_soluble'} ne 'all' ) {
					if ( (is_valid_positive_integer_from_zero($global->{'num_soluble'})) == 0 ) {
						$global->{'got_error'} = 1;
						$global->{'err_msg'} .= "The requested number of soluble sequences to include : " . $global->{'num_soluble'} . " is not a positive integer number.<br>";
					}
				}
			}
		}
	}

	if (($global->{'got_error'} == 0) and ($global->{'manual_list_seqs'} eq '')) {

		my @protein_categories = ( 'alpha_polytopic_bitopic', 'beta_barrel' );
		my @protein_categories_name = ( 'transmembrane helix (TMH)', 'beta-barrel', 'soluble protein' );
		my $i = 0;
		for my $this_protein_category ( @protein_categories ) {

			if ($global->{'kingdom'}->{$this_protein_category} ne '') {

				if (($global->{'similarity_percent'}->{$this_protein_category} ne '') and ($global->{'water_needle'}->{$this_protein_category} ne '') and ($global->{'similarity_identity'}->{$this_protein_category} ne '')) {
					my $input_params_have_valid_values = 0;
					my @valid_similarity_percent = (20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,101,'20_2','25_2','30_2','35_2','40_2','45_2','50_2','55_2','60_2','65_2','70_2','75_2','80_2','85_2','90_2','95_2','100_2');
					if (	(($global->{'water_needle'}->{$this_protein_category} eq 'water') or ($global->{'water_needle'}->{$this_protein_category} eq 'needle')) and
						(($global->{'similarity_identity'}->{$this_protein_category} eq 'similarity') or ($global->{'similarity_identity'}->{$this_protein_category} eq 'identity'))	) {
						foreach my $this_similarity_percent (@valid_similarity_percent) {
							if ($this_similarity_percent eq $global->{'similarity_percent'}->{$this_protein_category}) {
								$input_params_have_valid_values = 1;
							}
						}
					}
					my $found_this_data_file_in_list = 0;
					if ($input_params_have_valid_values == 1) {
						$global->{'nonhomologous_list'}->{$this_protein_category} = "data/$this_protein_category.pdb_seq." . $global->{'water_needle'}->{$this_protein_category} . '.' . $global->{'similarity_identity'}->{$this_protein_category} . '_scores.' . $global->{'similarity_percent'}->{$this_protein_category} . '.binary_matrix.nonhomologous_list';
						if ( index( $global->{'similarity_percent'}->{$this_protein_category}, '_2' ) > -1 ) {
							$global->{'similarity_percent'}->{$this_protein_category} =~ s/\_2//g;
							$global->{'nonhomologous_list'}->{$this_protein_category} = "data/$this_protein_category.pdb_seq.sort2." . $global->{'water_needle'}->{$this_protein_category} . '.' . $global->{'similarity_identity'}->{$this_protein_category} . '_scores.' . $global->{'similarity_percent'}->{$this_protein_category} . '.binary_matrix.nonhomologous_list';
						}
						my $file_having_list_of_data_files = $data->{'file_name'}->{$this_protein_category}->{'list_of_homology_files'};
						open LISTFILE, $file_having_list_of_data_files or die $!;
						my @input_lines = <LISTFILE>;
						foreach my $input_line (@input_lines) {
							chomp $input_line;
							if ($input_line ne '') {
								$input_line = trim($input_line);
								$input_line = "data/$input_line";
								if ($input_line eq $global->{'nonhomologous_list'}->{$this_protein_category}) {
									$found_this_data_file_in_list = 1;
								}
							}
						}
					}
					if ($found_this_data_file_in_list == 0) {
						$global->{'similarity_percent'}->{$this_protein_category} = '';
						$global->{'water_needle'}->{$this_protein_category} = '';
						$global->{'similarity_identity'}->{$this_protein_category} = '';
						$global->{'err_msg'} .= "The " . $protein_categories_name[$i] . " sequences chosen for benchmarking are not valid.<br>This is probably a system error.<br>\n";
						$global->{'input_benchmark_data_not_valid'} = 1;
						$global->{'got_error'} = 1;
					}
				}
			}
			$i++;
		}
	}

	if ($global->{'got_error'} == 0) {
		if ($global->{'run'} eq 'benchmark') {
			if (defined($params->{'minimum_overlap'})) {
				$params->{'minimum_overlap'} = trim($params->{'minimum_overlap'});
				if ( is_valid_positive_integer_above_zero($params->{'minimum_overlap'}) ) {
					$global->{'minimum_overlap'} = trim( $params->{'minimum_overlap'} );
				} else {
					$global->{'info_msg'} .= "Requested minimum overlap not a valid number and was adjusted to " . $data->{'minimum_overlap'} . ".<br>";
				}
			}
			if (defined($params->{'minimum_helix_length'})) {
				$params->{'minimum_helix_length'} = trim($params->{'minimum_helix_length'});
				if ( is_valid_positive_integer_above_zero($params->{'minimum_helix_length'}) ) {
					$global->{'minimum_helix_length'} = trim( $params->{'minimum_helix_length'} );
				} else {
					$global->{'info_msg'} .= "Requested minimum helix length not a valid number and was adjusted to " . $data->{'minimum_helix_length'} . ".<br>";
				}
			}
			if (defined($params->{'max_dist_approx_ends'})) {
				$params->{'max_dist_approx_ends'} = trim($params->{'max_dist_approx_ends'});
				if ( is_valid_positive_integer_from_zero($params->{'max_dist_approx_ends'}) ) {
					$global->{'max_dist_approx_ends'} = trim( $params->{'max_dist_approx_ends'} );
				} else {
					$global->{'info_msg'} .= "Requested maximum distance in residues from helix ends not a valid number and was adjusted to " . $data->{'max_dist_approx_ends'} . ".<br>";
				}
			}

			my @list_of_benchmark_criteria = @{$data->{'list_of_benchmark_criteria'}};
			my @list_of_benchmark_criteria_name = @{$data->{'list_of_benchmark_criteria_name'}};

			if (defined($params->{'sort_by'})) {
				$params->{'sort_by'} = trim($params->{'sort_by'});
				my $i = 0;
				foreach my $this_valid_sort_by (@list_of_benchmark_criteria) {
					if ($this_valid_sort_by eq $params->{'sort_by'}) {
						$global->{'sort_by'} = $params->{'sort_by'};
						$global->{'sort_by_name'} = $list_of_benchmark_criteria_name[$i];
						if ($global->{'sort_by'} eq 'none') {
							$global->{'sort_by'} = $data->{'sort_by'};
							$global->{'sort_by_name'} = $params->{'sort_by_name'};
						}
					}
					$i++;
				}
			}
			if (defined($params->{'graph_by'})) {
				$params->{'graph_by'} = trim($params->{'graph_by'});
				my $i = 0;
				foreach my $this_valid_sort_by (@list_of_benchmark_criteria) {
					if ($this_valid_sort_by eq $params->{'graph_by'}) {
						$global->{'graph_by'} = $params->{'graph_by'};
						if (($global->{'graph_by'} eq 'none') or ($global->{'graph_by'} eq 'alphabet')) {
							$global->{'graph_by'} = '';
						}
					}
					$i++;
				}
			}
			if (defined($params->{'graph_by_2'})) {
				$params->{'graph_by_2'} = trim($params->{'graph_by_2'});
				my $i = 0;
				foreach my $this_valid_sort_by (@list_of_benchmark_criteria) {
					if ($this_valid_sort_by eq $params->{'graph_by_2'}) {
						$global->{'graph_by_2'} = $params->{'graph_by_2'};
						if (($global->{'graph_by_2'} eq 'none') or ($global->{'graph_by_2'} eq 'alphabet')) {
							$global->{'graph_by_2'} = '';
						}
					}
					$i++;
				}
				if (($global->{'graph_by'} eq '') and ($global->{'graph_by_2'} ne '')) {
					$global->{'graph_by'} = $global->{'graph_by_2'};
					$global->{'graph_by_2'} = '';
				}
			}
			if (defined($params->{'graph_by_3'})) {
				$params->{'graph_by_3'} = trim($params->{'graph_by_3'});
				my $i = 0;
				foreach my $this_valid_sort_by (@list_of_benchmark_criteria) {
					if ($this_valid_sort_by eq $params->{'graph_by_3'}) {
						$global->{'graph_by_3'} = $params->{'graph_by_3'};
						if (($global->{'graph_by_3'} eq 'none') or ($global->{'graph_by_3'} eq 'alphabet')) {
							$global->{'graph_by_3'} = '';
						}
					}
					$i++;
				}
				if (($global->{'graph_by_2'} eq '') and ($global->{'graph_by_3'} ne '')) {
					$global->{'graph_by_2'} = $global->{'graph_by_3'};
					$global->{'graph_by_3'} = '';
				}
			}
			if (defined($params->{'graph_title'})) {
				$global->{'graph_title'} = trim($params->{'graph_title'});
			}
		}
	}

	if ($global->{'got_error'} == 0) {
		if (($global->{'run'} eq 'display') or ($global->{'run'} eq 'display_date')) {
			if (defined($params->{'display_width'})) {
				$params->{'display_width'} = trim($params->{'display_width'});
				if ( is_valid_positive_integer_from_zero($params->{'display_width'}) ) {
					$global->{'display_width'} = trim( $params->{'display_width'} );
				} else {
					$global->{'info_msg'} .= "Requested display width not a valid number and was adjusted to " . $data->{'display_width'} . ".<br>";
				}
			}
		}
	}

	# if the user uploaded a file, then read and validate the contents

	if ($global->{'got_error'} == 0) {

		my $upload_file = $q->param('upload_file');
		my $upload_filehandle = $q->upload('upload_file');
		if ($upload_file ne '') {

			$global->{'user_uploaded_file'} = 1;

			if ($upload_filehandle eq '') {
				$global->{'got_error'} = 1;
				$global->{'err_msg'} .= "The upload didn't work. Probably your file is too big - bigger than 5MB.<br>";
			} else {

				# read the uploaded file

				my $user_upload_file_data = '';
				while ( <$upload_filehandle> ) {
					my $input_line = $_;
					chomp($input_line);
					$user_upload_file_data .= $input_line . "\n";
				}

				read_and_validate_user_predictions( $user_upload_file_data );

			}
		}
	}

	if ($global->{'got_error'} == 0) {

	}

	# if the user entered predictions in the text box, then validate and save them,
	# but only if predictions were not already uploaded in a file.

	if ($global->{'got_error'} == 0) {

		my $user_predictions = $q->param('user_predictions');
		$user_predictions = trim($user_predictions);
		if ($user_predictions ne '') {
			$global->{'user_entered_predictions_in_text_box'} = 1;
			if ($global->{'user_uploaded_file'} == 1) {
				$global->{'user_file_and_textbox_predictions_present'} = 1;
			} else {

				read_and_validate_user_predictions( $user_predictions );

			}
		}
	}

	if ($global->{'got_error'} == 0) {

		my %output_display_method;
		$output->{'display_method'} = \%output_display_method;

		$output->{'display_method'}->{'opm_tm_segments'} = 0;
		$output->{'display_method'}->{'pdbtm_segments'} = 0;
		$output->{'display_method'}->{'opm_tmh_stride_segments'} = 0;
		$output->{'display_method'}->{'dssp'} = 0;
		$output->{'display_method'}->{'stride'} = 0;
		$output->{'display_method'}->{'opm_adj_tm_segments'} = 0;
		$output->{'display_method'}->{'opm_adj_tm_segments_topology'} = 0;
		$output->{'display_method'}->{'opm_tm_segments_topology'} = 0;
		$output->{'display_method'}->{'pdbtm_segments_topology'} = 0;

		if (defined($params->{'display_opm_adj'})) {
			$output->{'display_method'}->{'opm_adj_tm_segments'} = 1;
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
		}
		if (defined($params->{'display_opm_adj_topology'})) {
			$output->{'display_method'}->{'opm_adj_tm_segments_topology'} = 1;
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
		}
		if (defined($params->{'display_opm'})) {
			$output->{'display_method'}->{'opm_tm_segments'} = 1;
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
		}
		if (defined($params->{'display_opm_topology'})) {
			$output->{'display_method'}->{'opm_tm_segments_topology'} = 1;
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
		}
		if (defined($params->{'display_pdbtm'})) {
			$output->{'display_method'}->{'pdbtm_segments'} = 1;
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
		}
		if (defined($params->{'display_pdbtm_topology'})) {
			$output->{'display_method'}->{'pdbtm_segments_topology'} = 1;
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
		}
		if (defined($params->{'display_dssp'})) {
			$output->{'display_method'}->{'dssp'} = 1;
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
		}
		if (defined($params->{'display_stride'})) {
			$output->{'display_method'}->{'stride'} = 1;
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
		}
		if (defined($params->{'display_intermediate'})) {
			$global->{'show_intermediate'} = 1;
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 7;
		}
		if (defined($params->{'display_topology_files'})) {
			$global->{'show_topology_files'} = 1;
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 3;
		}
		if (defined($params->{'dont_show_sequences'})) {
			$global->{'dont_show_sequences'} = 1;
		}
		#----------
		$global->{'dastmfilter'} = 0;
		if (defined($params->{'dastmfilter'})) {
			$global->{'dastmfilter'} = trim( $params->{'dastmfilter'} );
		}
		if ($global->{'dastmfilter'} == 1) {
			push( @{$output->{'prediction_method'}}, 'dastmfilter' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'das2002'} = 0;
		if (defined($params->{'das2002'})) {
			$global->{'das2002'} = trim( $params->{'das2002'} );
		}
		if ($global->{'das2002'} == 1) {
			push( @{$output->{'prediction_method'}}, 'das2002' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'das1997_loose'} = 0;
		if (defined($params->{'das1997_loose'})) {
			$global->{'das1997_loose'} = trim( $params->{'das1997_loose'} );
		}
		if ($global->{'das1997_loose'} == 1) {
			push( @{$output->{'prediction_method'}}, 'das1997_loose' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'das1997_strict'} = 0;
		if (defined($params->{'das1997_strict'})) {
			$global->{'das1997_strict'} = trim( $params->{'das1997_strict'} );
		}
		if ($global->{'das1997_strict'} == 1) {
			push( @{$output->{'prediction_method'}}, 'das1997_strict' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'deltag'} = 0;
		if (defined($params->{'deltag'})) {
			$global->{'deltag'} = trim( $params->{'deltag'} );
		}
		if ($global->{'deltag'} == 1) {
			push( @{$output->{'prediction_method'}}, 'deltag' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'ensemble'} = 0;
		if (defined($params->{'ensemble'})) {
			$global->{'ensemble'} = trim( $params->{'ensemble'} );
		}
		if ($global->{'ensemble'} == 1) {
			push( @{$output->{'prediction_method'}}, 'ensemble' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'pepinfo_eisenberg_7_10'} = 0;
		if (defined($params->{'pepinfo_eisenberg_7_10'})) {
			$global->{'pepinfo_eisenberg_7_10'} = trim( $params->{'pepinfo_eisenberg_7_10'} );
		}
		if ($global->{'pepinfo_eisenberg_7_10'} == 1) {
			push( @{$output->{'prediction_method'}}, 'pepinfo_eisenberg_7_10' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'pepinfo_eisenberg_11_10'} = 0;
		if (defined($params->{'pepinfo_eisenberg_11_10'})) {
			$global->{'pepinfo_eisenberg_11_10'} = trim( $params->{'pepinfo_eisenberg_11_10'} );
		}
		if ($global->{'pepinfo_eisenberg_11_10'} == 1) {
			push( @{$output->{'prediction_method'}}, 'pepinfo_eisenberg_11_10' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'pepinfo_eisenberg_19_10'} = 0;
		if (defined($params->{'pepinfo_eisenberg_19_10'})) {
			$global->{'pepinfo_eisenberg_19_10'} = trim( $params->{'pepinfo_eisenberg_19_10'} );
		}
		if ($global->{'pepinfo_eisenberg_19_10'} == 1) {
			push( @{$output->{'prediction_method'}}, 'pepinfo_eisenberg_19_10' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'hmmtm'} = 0;
		if (defined($params->{'hmmtm'})) {
			$global->{'hmmtm'} = trim( $params->{'hmmtm'} );
		}
		if ($global->{'hmmtm'} == 1) {
			push( @{$output->{'prediction_method'}}, 'hmmtm' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'hmmtop2'} = 0;
		if (defined($params->{'hmmtop2'})) {
			$global->{'hmmtop2'} = trim( $params->{'hmmtop2'} );
		}
		if ($global->{'hmmtop2'} == 1) {
			push( @{$output->{'prediction_method'}}, 'hmmtop2' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'hmmtop_topcons_single'} = 0;
		if (defined($params->{'hmmtop_topcons_single'})) {
			$global->{'hmmtop_topcons_single'} = trim( $params->{'hmmtop_topcons_single'} );
		}
		if ($global->{'hmmtop_topcons_single'} == 1) {
			push( @{$output->{'prediction_method'}}, 'hmmtop_topcons_single' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'pepinfo_kyte_7_10'} = 0;
		if (defined($params->{'pepinfo_kyte_7_10'})) {
			$global->{'pepinfo_kyte_7_10'} = trim( $params->{'pepinfo_kyte_7_10'} );
		}
		if ($global->{'pepinfo_kyte_7_10'} == 1) {
			push( @{$output->{'prediction_method'}}, 'pepinfo_kyte_7_10' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'pepinfo_kyte_11_10'} = 0;
		if (defined($params->{'pepinfo_kyte_11_10'})) {
			$global->{'pepinfo_kyte_11_10'} = trim( $params->{'pepinfo_kyte_11_10'} );
		}
		if ($global->{'pepinfo_kyte_11_10'} == 1) {
			push( @{$output->{'prediction_method'}}, 'pepinfo_kyte_11_10' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'pepinfo_kyte_19_10'} = 0;
		if (defined($params->{'pepinfo_kyte_19_10'})) {
			$global->{'pepinfo_kyte_19_10'} = trim( $params->{'pepinfo_kyte_19_10'} );
		}
		if ($global->{'pepinfo_kyte_19_10'} == 1) {
			push( @{$output->{'prediction_method'}}, 'pepinfo_kyte_19_10' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'membrain'} = 0;
		if (defined($params->{'membrain'})) {
			$global->{'membrain'} = trim( $params->{'membrain'} );
		}
		if ($global->{'membrain'} == 1) {
			push( @{$output->{'prediction_method'}}, 'membrain' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'memsatsvm'} = 0;
		if (defined($params->{'memsatsvm'})) {
			$global->{'memsatsvm'} = trim( $params->{'memsatsvm'} );
		}
		if ($global->{'memsatsvm'} == 1) {
			push( @{$output->{'prediction_method'}}, 'memsatsvm' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'memsat3'} = 0;
		if (defined($params->{'memsat3'})) {
			$global->{'memsat3'} = trim( $params->{'memsat3'} );
		}
		if ($global->{'memsat3'} == 1) {
			push( @{$output->{'prediction_method'}}, 'memsat3' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'memsat_topcons_single'} = 0;
		if (defined($params->{'memsat_topcons_single'})) {
			$global->{'memsat_topcons_single'} = trim( $params->{'memsat_topcons_single'} );
		}
		if ($global->{'memsat_topcons_single'} == 1) {
			push( @{$output->{'prediction_method'}}, 'memsat_topcons_single' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'minnou'} = 0;
		if (defined($params->{'minnou'})) {
			$global->{'minnou'} = trim( $params->{'minnou'} );
		}
		if ($global->{'minnou'} == 1) {
			push( @{$output->{'prediction_method'}}, 'minnou' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'octopus'} = 0;
		if (defined($params->{'octopus'})) {
			$global->{'octopus'} = trim( $params->{'octopus'} );
		}
		if ($global->{'octopus'} == 1) {
			push( @{$output->{'prediction_method'}}, 'octopus' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'octopus_topcons'} = 0;
		if (defined($params->{'octopus_topcons'})) {
			$global->{'octopus_topcons'} = trim( $params->{'octopus_topcons'} );
		}
		if ($global->{'octopus_topcons'} == 1) {
			push( @{$output->{'prediction_method'}}, 'octopus_topcons' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'pepinfo_ohm_7_10'} = 0;
		if (defined($params->{'pepinfo_ohm_7_10'})) {
			$global->{'pepinfo_ohm_7_10'} = trim( $params->{'pepinfo_ohm_7_10'} );
		}
		if ($global->{'pepinfo_ohm_7_10'} == 1) {
			push( @{$output->{'prediction_method'}}, 'pepinfo_ohm_7_10' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'pepinfo_ohm_11_10'} = 0;
		if (defined($params->{'pepinfo_ohm_11_10'})) {
			$global->{'pepinfo_ohm_11_10'} = trim( $params->{'pepinfo_ohm_11_10'} );
		}
		if ($global->{'pepinfo_ohm_11_10'} == 1) {
			push( @{$output->{'prediction_method'}}, 'pepinfo_ohm_11_10' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'pepinfo_ohm_19_10'} = 0;
		if (defined($params->{'pepinfo_ohm_19_10'})) {
			$global->{'pepinfo_ohm_19_10'} = trim( $params->{'pepinfo_ohm_19_10'} );
		}
		if ($global->{'pepinfo_ohm_19_10'} == 1) {
			push( @{$output->{'prediction_method'}}, 'pepinfo_ohm_19_10' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'phdhtm_pbil'} = 0;
		if (defined($params->{'phdhtm_pbil'})) {
			$global->{'phdhtm_pbil'} = trim( $params->{'phdhtm_pbil'} );
		}
		if ($global->{'phdhtm_pbil'} == 1) {
			push( @{$output->{'prediction_method'}}, 'phdhtm_pbil' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'phdthtm_pbil'} = 0;
		if (defined($params->{'phdthtm_pbil'})) {
			$global->{'phdthtm_pbil'} = trim( $params->{'phdthtm_pbil'} );
		}
		if ($global->{'phdthtm_pbil'} == 1) {
			push( @{$output->{'prediction_method'}}, 'phdthtm_pbil' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'philius'} = 0;
		if (defined($params->{'philius'})) {
			$global->{'philius'} = trim( $params->{'philius'} );
		}
		if ($global->{'philius'} == 1) {
			push( @{$output->{'prediction_method'}}, 'philius' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'phobius'} = 0;
		if (defined($params->{'phobius'})) {
			$global->{'phobius'} = trim( $params->{'phobius'} );
		}
		if ($global->{'phobius'} == 1) {
			push( @{$output->{'prediction_method'}}, 'phobius' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'polyphobius'} = 0;
		if (defined($params->{'polyphobius'})) {
			$global->{'polyphobius'} = trim( $params->{'polyphobius'} );
		}
		if ($global->{'polyphobius'} == 1) {
			push( @{$output->{'prediction_method'}}, 'polyphobius' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'predtmr'} = 0;
		if (defined($params->{'predtmr'})) {
			$global->{'predtmr'} = trim( $params->{'predtmr'} );
		}
		if ($global->{'predtmr'} == 1) {
			push( @{$output->{'prediction_method'}}, 'predtmr' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'pro_topcons'} = 0;
		if (defined($params->{'pro_topcons'})) {
			$global->{'pro_topcons'} = trim( $params->{'pro_topcons'} );
		}
		if ($global->{'pro_topcons'} == 1) {
			push( @{$output->{'prediction_method'}}, 'pro_topcons' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'prodiv_topcons'} = 0;
		if (defined($params->{'prodiv_topcons'})) {
			$global->{'prodiv_topcons'} = trim( $params->{'prodiv_topcons'} );
		}
		if ($global->{'prodiv_topcons'} == 1) {
			push( @{$output->{'prediction_method'}}, 'prodiv_topcons' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'scampi'} = 0;
		if (defined($params->{'scampi'})) {
			$global->{'scampi'} = trim( $params->{'scampi'} );
		}
		if ($global->{'scampi'} == 1) {
			push( @{$output->{'prediction_method'}}, 'scampi' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'scampi_msa_topcons'} = 0;
		if (defined($params->{'scampi_msa_topcons'})) {
			$global->{'scampi_msa_topcons'} = trim( $params->{'scampi_msa_topcons'} );
		}
		if ($global->{'scampi_msa_topcons'} == 1) {
			push( @{$output->{'prediction_method'}}, 'scampi_msa_topcons' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'scampi_seq_topcons'} = 0;
		if (defined($params->{'scampi_seq_topcons'})) {
			$global->{'scampi_seq_topcons'} = trim( $params->{'scampi_seq_topcons'} );
		}
		if ($global->{'scampi_seq_topcons'} == 1) {
			push( @{$output->{'prediction_method'}}, 'scampi_seq_topcons' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'scampi_seq_topcons_single'} = 0;
		if (defined($params->{'scampi_seq_topcons_single'})) {
			$global->{'scampi_seq_topcons_single'} = trim( $params->{'scampi_seq_topcons_single'} );
		}
		if ($global->{'scampi_seq_topcons_single'} == 1) {
			push( @{$output->{'prediction_method'}}, 'scampi_seq_topcons_single' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'sosui'} = 0;
		if (defined($params->{'sosui'})) {
			$global->{'sosui'} = trim( $params->{'sosui'} );
		}
		if ($global->{'sosui'} == 1) {
			push( @{$output->{'prediction_method'}}, 'sosui' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'split4'} = 0;
		if (defined($params->{'split4'})) {
			$global->{'split4'} = trim( $params->{'split4'} );
		}
		if ($global->{'split4'} == 1) {
			push( @{$output->{'prediction_method'}}, 'split4' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'svmtm'} = 0;
		if (defined($params->{'svmtm'})) {
			$global->{'svmtm'} = trim( $params->{'svmtm'} );
		}
		if ($global->{'svmtm'} == 1) {
			push( @{$output->{'prediction_method'}}, 'svmtm' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'svmtop'} = 0;
		if (defined($params->{'svmtop'})) {
			$global->{'svmtop'} = trim( $params->{'svmtop'} );
		}
		if ($global->{'svmtop'} == 1) {
			push( @{$output->{'prediction_method'}}, 'svmtop' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'tmap'} = 0;
		if (defined($params->{'tmap'})) {
			$global->{'tmap'} = trim( $params->{'tmap'} );
		}
		if ($global->{'tmap'} == 1) {
			push( @{$output->{'prediction_method'}}, 'tmap' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'tmhmm2'} = 0;
		if (defined($params->{'tmhmm2'})) {
			$global->{'tmhmm2'} = trim( $params->{'tmhmm2'} );
		}
		if ($global->{'tmhmm2'} == 1) {
			push( @{$output->{'prediction_method'}}, 'tmhmm2' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'stmhmm_topcons_single'} = 0;
		if (defined($params->{'stmhmm_topcons_single'})) {
			$global->{'stmhmm_topcons_single'} = trim( $params->{'stmhmm_topcons_single'} );
		}
		if ($global->{'stmhmm_topcons_single'} == 1) {
			push( @{$output->{'prediction_method'}}, 'stmhmm_topcons_single' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'tmloop'} = 0;
		if (defined($params->{'tmloop'})) {
			$global->{'tmloop'} = trim( $params->{'tmloop'} );
		}
		if ($global->{'tmloop'} == 1) {
			push( @{$output->{'prediction_method'}}, 'tmloop' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'tmmod'} = 0;
		if (defined($params->{'tmmod'})) {
			$global->{'tmmod'} = trim( $params->{'tmmod'} );
		}
		if ($global->{'tmmod'} == 1) {
			push( @{$output->{'prediction_method'}}, 'tmmod' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'tmpred'} = 0;
		if (defined($params->{'tmpred'})) {
			$global->{'tmpred'} = trim( $params->{'tmpred'} );
		}
		if ($global->{'tmpred'} == 1) {
			push( @{$output->{'prediction_method'}}, 'tmpred' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'topcons'} = 0;
		if (defined($params->{'topcons'})) {
			$global->{'topcons'} = trim( $params->{'topcons'} );
		}
		if ($global->{'topcons'} == 1) {
			push( @{$output->{'prediction_method'}}, 'topcons' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'topcons_single'} = 0;
		if (defined($params->{'topcons_single'})) {
			$global->{'topcons_single'} = trim( $params->{'topcons_single'} );
		}
		if ($global->{'topcons_single'} == 1) {
			push( @{$output->{'prediction_method'}}, 'topcons_single' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'toppred2'} = 0;
		if (defined($params->{'toppred2'})) {
			$global->{'toppred2'} = trim( $params->{'toppred2'} );
		}
		if ($global->{'toppred2'} == 1) {
			push( @{$output->{'prediction_method'}}, 'toppred2' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'valpred'} = 0;
		if (defined($params->{'valpred'})) {
			$global->{'valpred'} = trim( $params->{'valpred'} );
		}
		if ($global->{'valpred'} == 1) {
			push( @{$output->{'prediction_method'}}, 'valpred' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'valpred2'} = 0;
		if (defined($params->{'valpred2'})) {
			$global->{'valpred2'} = trim( $params->{'valpred2'} );
		}
		if ($global->{'valpred2'} == 1) {
			push( @{$output->{'prediction_method'}}, 'valpred2' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
		$global->{'wavetm'} = 0;
		if (defined($params->{'wavetm'})) {
			$global->{'wavetm'} = trim( $params->{'wavetm'} );
		}
		if ($global->{'wavetm'} == 1) {
			push( @{$output->{'prediction_method'}}, 'wavetm' );
			$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
			$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;
		}
		#----------
	}

	if ($global->{'got_error'} == 0) {
		if ($global->{'run'} eq 'benchmark') {
			if ($global->{'num_output_prediction_methods'} == 0) {
				$global->{'err_msg'} .= "No prediction methods were chosen.<br>Please choose at least one method for benchmarking.<br>\n";
				$global->{'got_error'} = 1;	
			}
		}
	}
}

#===================================================================================================
#	read in and validate the user predictions, 
#	which have come from uploaded file or were entered in html form text box.
#===================================================================================================

sub read_and_validate_user_predictions {

	my $user_data = shift;
	my @user_data_lines = split( /\n|\r/, $user_data );

	my %user_data_hash;
	$data->{'user_data'} = \%user_data_hash;

	my $seqid = '';
	my $this_segments_sequence = '';
	my $got_prev_segments = 0;

	$global->{'use_the_user_predictions'} = 1;
	my $all_sequence_lines_are_valid = 1;
	my $format_of_lines_are_all_valid = 1;
	my $found_a_zero = 0;
	my $prev_line_type = 'seq'; # 'id', 'seq'
	my $curr_line_type = '';

	foreach my $input_line (@user_data_lines) {
		chomp($input_line);
		$input_line = trim($input_line);
		if (($input_line ne '') and ($global->{'use_the_user_predictions'} == 1)) {

			if (substr($input_line,0,1) eq '>') {

				# save the previous id
				if ($got_prev_segments == 1) {
					if ($this_segments_sequence eq '') {
						$this_segments_sequence = '-';
					}
					$data->{'user_data'}->{$seqid} = $this_segments_sequence;
					$got_prev_segments = 0;
				}

				$seqid = substr($input_line, 1);
				$this_segments_sequence = '';
				$got_prev_segments = 1;
				$prev_line_type = 'id';

			} else {

				if ($prev_line_type eq 'seq') {
					if ($global->{'use_the_user_predictions'} == 1) {
						$format_of_lines_are_all_valid = 0;
					}

				} else {

					my $this_segments_line = $input_line;
					$this_segments_sequence = '';

					my $ix = index($this_segments_line, ',');

					if ((index($this_segments_line, ',') > -1) || (index($this_segments_line, '-') > -1)) {

						# 11-30,53-71
						my @bits = split( /\,/, $this_segments_line );
						my $prev_end = 0;
						foreach my $pair (@bits) {
							if ($pair ne '') {
								if ($all_sequence_lines_are_valid == 1) {
									$pair =~ /(-?\d+)-(-?\d+)/;
									if ( (!defined($1)) || (!defined($2)) ) {
										$all_sequence_lines_are_valid = 0;
									} else {
										my $start = $1;
										my $end = $2;
										if (($start !~ /^\d+$/) || ($end !~ /^\d+$/)) {
											$all_sequence_lines_are_valid = 0;
										} else {
											if (($start == 0) || ($end == 0)) {
												$found_a_zero = 1;
											} else {
												for ( my $i = ($prev_end+1); $i <= ($start-1); $i++ ) {
													$this_segments_sequence .= '_';
												}
												for ( my $i = $start; $i <= $end; $i++ ) {
													$this_segments_sequence .= 'M';
												}
											}
										}
										$prev_end = $end;
									}
								}
							}
						}

					} else {

						# ooooooooooMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMooooo
						$this_segments_sequence = $this_segments_line;
					}
				}
				$prev_line_type = 'seq';
			}
		}
	}
	# save the last id
	if ($got_prev_segments == 1) {
		if ($this_segments_sequence eq '') {
			$this_segments_sequence = '-';
		}
		$data->{'user_data'}->{$seqid} = $this_segments_sequence;
	}

	if (($got_prev_segments == 1) and ($global->{'use_the_user_predictions'} == 1)) {
		if ($this_segments_sequence eq '') {
			$this_segments_sequence = '-';
		}
		$data->{'user_data'}->{$seqid} = $this_segments_sequence;
		$got_prev_segments = 0;
	}

	if ($format_of_lines_are_all_valid == 0) {
		$global->{'got_error'} = 1;
		$global->{'err_msg'} .= "User-entered prediction lines do not seem to be in the correct format.<br>";
		$global->{'use_the_user_predictions'} = 0;
	}
	if ($all_sequence_lines_are_valid == 0) {
		$global->{'got_error'} = 1;
		$global->{'err_msg'} .= "User-entered prediction numbers do not seem to be in the correct format.<br>";
		$global->{'use_the_user_predictions'} = 0;
	}
	if ($found_a_zero == 1) {
		$global->{'got_error'} = 1;
		$global->{'err_msg'} .= "Make sure your segments start at position 1, not position 0.<br>";
		$global->{'use_the_user_predictions'} = 0;
	}

	foreach my $seqid ( keys %{$data->{'user_data'}} ) {
		if ( exists($data->{'user_data'}->{$seqid}) ) {
			$output->{'num_user_seqs'} = $output->{'num_user_seqs'} + 1;
			my $segments = $data->{'user_data'}->{$seqid};
			if ($segments ne '-') {
				my @bits = split( /\,/, $segments );
				foreach my $this_segment (@bits) {
					$output->{'num_user_membrane_helices'} = $output->{'num_user_membrane_helices'} + 1;
				}
			}
		}
	}

	if ($global->{'use_the_user_predictions'} == 1) {

		push( @{$output->{'prediction_method'}}, 'you' );
		$global->{'num_output_display_methods'} = $global->{'num_output_display_methods'} + 1;
		$global->{'num_output_prediction_methods'} = $global->{'num_output_prediction_methods'} + 1;

		my $has_topology = 0;
		foreach my $this_user_data ( values %{$data->{'user_data'}} ) {
			if ($this_user_data =~ m/i/) {
				$has_topology = 1;
			}
			if ($this_user_data =~ m/o/) {
				$has_topology = 1;
			}
		}
		if ($has_topology == 1) {
			$global->{'input_benchmark_data_has_topology'} = 1;
			push( @{$data->{'list_of_topology_prediction_methods'}}, 'you' );
		}
	}
}

#===================================================================================================
#	now that we have the list of sequences for benchmarking, 
#	see if this matches the benchmark predictions that the user entered.
#===================================================================================================

sub validate_user_predictions_list {

	if ($global->{'use_the_user_predictions'} == 1) {

		foreach my $seqid ( keys %{$output->{'id'}->{'all'}} ) {
			if (!defined( $data->{'user_data'}->{$seqid} )) {
				$global->{'user_predictions_missing_defined_sequences'} = 1;
			}
		}
		foreach my $seqid ( keys %{$data->{'user_data'}} ) {
			my $found_seqid = 0;
			foreach my $output_seqid ( keys %{$output->{'id'}->{'all'}} ) {
				if ($output_seqid eq $seqid) {
					$found_seqid = 1;
				}
			}
			if ($found_seqid == 0) {
				$global->{'user_predictions_have_extra_sequences'} = 1;
			}
		}
	}
}

#===================================================================================================
#	now that we have the sequences, validate that user entered predictions are not too long
#	also pad out the length of the sequence if necessary
#	(if user predictions are in the format 11-30,53-71, will need to pad out the length for display now that we have read the sequence)
#===================================================================================================

sub validate_user_prediction_sequence_lengths {

	my $user_data_hash = $data->{'user_data'};
	my %user_data = %$user_data_hash;

	foreach my $seqid (keys %user_data) {
		my $this_seq_length = length( $output->{'aaseq'}->{$seqid} );
		if ( exists($output->{'aaseq'}->{$seqid}) ) {
			my $this_user_length = length( $data->{'user_data'}->{$seqid} );
			if ($this_user_length > $this_seq_length) {
				if ($global->{'user_prediction_longer_than_sequence'} == 0) {
					$global->{'info_msg'} .= "A user-entered prediction is longer than the sequence length.<br>\n";
					$global->{'user_prediction_longer_than_sequence'} = 1;
				}
			} elsif ($this_user_length < $this_seq_length) {
				my $new_seq = $data->{'user_data'}->{$seqid};
				for ( my $i = ($this_user_length+1); $i <= $this_seq_length; $i++ ) {
					$new_seq .= '_';
				}
				$data->{'user_data'}->{$seqid} = $new_seq;
			}
		}
	}
}

#===================================================================================================
#	initialisations for this script
#===================================================================================================

sub init {

	$CGI::POST_MAX = 1024 * 5000; # limit file uploads to 5MB

	my %global_hash;
	my %data_hash;
	my %output_hash;
	$global = \%global_hash;
	$data = \%data_hash;
	$output = \%output_hash;

	$global->{'program_mode'} = '';
	$global->{'got_error'} = 0;
	$global->{'err_msg'} = '';
	$global->{'info_msg'} = '';
	$global->{'display_text_not_html'} = 0;
	$global->{'run'} = '';
	$global->{'user_predictions_missing_defined_sequences'} = 0;
	$global->{'user_predictions_have_extra_sequences'} = 0;
	$global->{'user_uploaded_file'} = 0;
	$global->{'user_entered_predictions_in_text_box'} = 0;
	$global->{'user_file_and_textbox_predictions_present'} = 0;
	$global->{'user_prediction_longer_than_sequence'} = 0;
	$global->{'use_the_user_predictions'} = 0;
	$global->{'minimum_overlap_not_valid'} = 0;
	$global->{'minimum_helix_length_not_valid'} = 0;
	$global->{'max_dist_approx_ends_not_valid'} = 0;
	$global->{'display_width_not_valid'} = 0;
	$global->{'input_benchmark_data_not_valid'} = 0;
	$global->{'input_benchmark_data_has_topology'} = 0;
	$global->{'benchmark_standard'} = '';
	$global->{'show_intermediate'} = 0;
	$global->{'show_topology_files'} = 0;
	$global->{'dont_show_sequences'} = 0;
	$global->{'manual_list_seqs'} = '';

	$data->{'minimum_overlap'} = 9;
	$data->{'minimum_helix_length'} = 5;
	$data->{'max_dist_approx_ends'} = 5;
	$data->{'display_width'} = 0;
	$data->{'sort_by'} = 'qok';
	$data->{'sort_by_name'} = 'Qok&#37;';

	my @list_of_benchmark_criteria =      ('qok',     'qhtm_obs',                   'qhtm_prd',                   'avhe_diff','qhe_obs',     'qhe_prd',     'qnhe_obs',          'qnhe_prd',          'q2',     'htm_mcc','q2t_obs',     'q2t_prd',     'q2n_obs',     'q2n_prd',     'qok3',     'nterm',      'q2_ioseg',     'qiom_obs',     'qiom_prd',     'qio_obs',     'qio_prd',     'q3',    'q2_iores',      'io_mcc','alphabet',   'none');
	my @list_of_benchmark_criteria_name = ('Qok&#37;','Qhtm &#37;obs (sensitivity)','Qhtm &#37;prd (specificity)','AvHb diff','QHb &#37;obs','QHb &#37;prd','Gauss QHb &#37;obs','Gauss QHb &#37;prd','Q2&#37;','htm MCC','Q2T &#37;obs','Q2T &#37;prd','Q2N &#37;obs','Q2N &#37;prd','Qok3&#37;','Nterm &#37;','ioSeg Q2&#37;','Qiom &#37;obs','Qiom &#37;prd','Qio &#37;obs','Qio &#37;prd','Q3&#37;','ioRes Q2&#37;','io MCC','Method Name','No sort');
	$data->{'list_of_benchmark_criteria'} = \@list_of_benchmark_criteria;
	$data->{'list_of_benchmark_criteria_name'} = \@list_of_benchmark_criteria_name;

	my @list_of_prediction_methods = ( 'dastmfilter', 'das2002', 'das1997_loose', 'das1997_strict', 'deltag', 'ensemble', 'hmmtm', 'hmmtop2', 'hmmtop_topcons_single', 'membrain', 'memsatsvm', 'memsat3', 'memsat_topcons_single', 'minnou', 'octopus', 'octopus_topcons', 'phdhtm_pbil', 'phdthtm_pbil', 'philius', 'phobius', 'polyphobius', 'predtmr', 'pro_topcons', 'prodiv_topcons', 'scampi', 'scampi_msa_topcons', 'scampi_seq_topcons', 'scampi_seq_topcons_single', 'sosui', 'split4', 'svmtm', 'svmtop', 'tmap', 'tmhmm2', 'stmhmm_topcons_single', 'tmloop', 'tmmod', 'tmpred', 'topcons', 'topcons_single', 'toppred2', 'valpred', 'valpred2', 'wavetm', 'pepinfo_eisenberg_7_10', 'pepinfo_eisenberg_11_10', 'pepinfo_eisenberg_19_10', 'pepinfo_kyte_7_10', 'pepinfo_kyte_11_10', 'pepinfo_kyte_19_10', 'pepinfo_ohm_7_10', 'pepinfo_ohm_11_10', 'pepinfo_ohm_19_10', 'you' );
	$data->{'list_of_prediction_methods'} = \@list_of_prediction_methods;
	my @list_of_topology_prediction_methods = ( 'ensemble', 'hmmtm', 'hmmtop2', 'hmmtop_topcons_single', 'memsatsvm', 'memsat3', 'memsat_topcons_single', 'octopus', 'octopus_topcons', 'phdthtm_pbil', 'philius', 'phobius', 'polyphobius', 'pro_topcons', 'prodiv_topcons', 'scampi', 'scampi_msa_topcons', 'scampi_seq_topcons', 'scampi_seq_topcons_single', 'svmtop', 'tmhmm2', 'stmhmm_topcons_single', 'tmmod', 'tmpred', 'topcons', 'topcons_single', 'toppred2' );
	$data->{'list_of_topology_prediction_methods'} = \@list_of_topology_prediction_methods;

	my %data_file_name;
	$data->{'file_name'} = \%data_file_name;
	my %data_file_name_alpha_polytopic_bitopic;
	$data->{'file_name'}->{'alpha_polytopic_bitopic'} = \%data_file_name_alpha_polytopic_bitopic;
	my %data_file_name_beta_barrel;
	$data->{'file_name'}->{'beta_barrel'} = \%data_file_name_beta_barrel;

	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'list_of_homology_files'} = 'data/alpha_polytopic_bitopic.nonhomologous_list.files';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'id'} = 'data/alpha_polytopic_bitopic.id';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'mpstruc_category'} = 'data/alpha_polytopic_bitopic.mpstruc_category';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'opm_category'} = 'data/alpha_polytopic_bitopic.opm_category';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'pdbe_name'} = 'data/alpha_polytopic_bitopic.pdbe_name';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'pdbe_taxonomy'} = 'data/alpha_polytopic_bitopic.pdbe_taxonomy';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'release_date'} = 'data/alpha_polytopic_bitopic.release_date';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'reentrant_helices'} = 'data/alpha_polytopic_bitopic.reentrant_helices';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'dssp'} = 'data/alpha_polytopic_bitopic.pdb_struc2D';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'opm_adj_tm_segments'} = 'data/alpha_polytopic_bitopic.opm_adj_tm_segments';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'opm_adj_tm_segments_topology'} = 'data/alpha_polytopic_bitopic.opm_adj_tm_segments_topology';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'opm_adj_tmh_tmbb_segments'} = 'data/alpha_polytopic_bitopic.opm_adj_tmh_tmbb_segments';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'opm_membrane_near'} = 'data/alpha_polytopic_bitopic.opm_membrane_near';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'opm_posn'} = 'data/alpha_polytopic_bitopic.opm_posn';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'opm_seq'} = 'data/alpha_polytopic_bitopic.opm_seq';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'opm_tm_segments'} = 'data/alpha_polytopic_bitopic.opm_tm_segments'; # includes membrane helices and beta-barrel as MMMMMM
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'opm_tm_segments_topology'} = 'data/alpha_polytopic_bitopic.opm_tm_segments_topology'; 
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'opm_tmh_tmbb_segments'} = 'data/alpha_polytopic_bitopic.opm_tmh_tmbb_segments'; # includes membrane helices and beta-barrel as HHHHHH or BBBBBB
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'opm_tm_subunits'} = 'data/alpha_polytopic_bitopic.opm_tm_subunits';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'pdb_seq'} = 'data/alpha_polytopic_bitopic.pdb_seq';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'pdbtm_segments'} = 'data/alpha_polytopic_bitopic.pdbtm_segments';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'pdbtm_segments_topology'} = 'data/alpha_polytopic_bitopic.pdbtm_segments_topology';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'stride'} = 'data/alpha_polytopic_bitopic.stride_struc2D';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'stride_seq'} = 'data/alpha_polytopic_bitopic.stride_seq';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'earliest_release_date'} = 'data/alpha_polytopic_bitopic.earliest_release_date';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'pdb_method_resolution'} = 'data/alpha_polytopic_bitopic.pdb_method_resolution';
	$data->{'file_name'}->{'beta_barrel'}->{'list_of_homology_files'} = 'data/beta_barrel.nonhomologous_list.files';
	$data->{'file_name'}->{'beta_barrel'}->{'id'} = 'data/beta_barrel.id';
	#$data->{'file_name'}->{'beta_barrel'}->{'mpstruc_category'} = 'data/beta_barrel.mpstruc_category';
	$data->{'file_name'}->{'beta_barrel'}->{'opm_category'} = 'data/beta_barrel.opm_category';
	$data->{'file_name'}->{'beta_barrel'}->{'pdbe_name'} = 'data/beta_barrel.pdbe_name';
	$data->{'file_name'}->{'beta_barrel'}->{'pdbe_taxonomy'} = 'data/beta_barrel.pdbe_taxonomy';
	$data->{'file_name'}->{'beta_barrel'}->{'release_date'} = 'data/beta_barrel.release_date';
	#$data->{'file_name'}->{'beta_barrel'}->{'reentrant_helices'} = 'data/beta_barrel.reentrant_helices';
	$data->{'file_name'}->{'beta_barrel'}->{'dssp'} = 'data/beta_barrel.pdb_struc2D';
	$data->{'file_name'}->{'beta_barrel'}->{'opm_adj_tm_segments'} = 'data/beta_barrel.opm_adj_tm_segments';
	$data->{'file_name'}->{'beta_barrel'}->{'opm_adj_tm_segments_topology'} = 'data/beta_barrel.opm_adj_tm_segments_topology';
	$data->{'file_name'}->{'beta_barrel'}->{'opm_adj_tmh_tmbb_segments'} = 'data/beta_barrel.opm_adj_tmh_tmbb_segments';
	$data->{'file_name'}->{'beta_barrel'}->{'opm_membrane_near'} = 'data/beta_barrel.opm_membrane_near';
	$data->{'file_name'}->{'beta_barrel'}->{'opm_posn'} = 'data/beta_barrel.opm_posn';
	$data->{'file_name'}->{'beta_barrel'}->{'opm_seq'} = 'data/beta_barrel.opm_seq';
	$data->{'file_name'}->{'beta_barrel'}->{'opm_tm_segments'} = 'data/beta_barrel.opm_tm_segments'; # includes membrane helices and beta-barrel as MMMMMM
	$data->{'file_name'}->{'beta_barrel'}->{'opm_tm_segments_topology'} = 'data/beta_barrel.opm_tm_segments_topology'; 
	$data->{'file_name'}->{'beta_barrel'}->{'opm_tmh_tmbb_segments'} = 'data/beta_barrel.opm_tmh_tmbb_segments'; # includes membrane helices and beta-barrel as HHHHHH or BBBBBB
	$data->{'file_name'}->{'beta_barrel'}->{'opm_tm_subunits'} = 'data/beta_barrel.opm_tm_subunits';
	$data->{'file_name'}->{'beta_barrel'}->{'pdb_seq'} = 'data/beta_barrel.pdb_seq';
	$data->{'file_name'}->{'beta_barrel'}->{'pdbtm_segments'} = 'data/beta_barrel.pdbtm_segments';
	$data->{'file_name'}->{'beta_barrel'}->{'pdbtm_segments_topology'} = 'data/beta_barrel.pdbtm_segments_topology';
	$data->{'file_name'}->{'beta_barrel'}->{'stride'} = 'data/beta_barrel.stride_struc2D';
	$data->{'file_name'}->{'beta_barrel'}->{'stride_seq'} = 'data/beta_barrel.stride_seq';
	$data->{'file_name'}->{'beta_barrel'}->{'pdb_method_resolution'} = 'data/beta_barrel.pdb_method_resolution';
	$data->{'file_name'}->{'soluble'}->{'id'} = 'data/soluble.id';
	$data->{'file_name'}->{'soluble'}->{'pdbe_name'} = 'data/soluble.pdbe_name';
	$data->{'file_name'}->{'soluble'}->{'pdbe_taxonomy'} = 'data/soluble.pdbe_taxonomy';
	$data->{'file_name'}->{'soluble'}->{'release_date'} = 'data/soluble.release_date';
	$data->{'file_name'}->{'soluble'}->{'dssp'} = 'data/soluble.pdb_struc2D';
	$data->{'file_name'}->{'soluble'}->{'pdb_seq'} = 'data/soluble.pdb_seq';
	$data->{'file_name'}->{'soluble'}->{'pdb_method_resolution'} = 'data/soluble.pdb_method_resolution';
	#----------
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'das1997_loose'} = 'data/alpha_polytopic_bitopic.das1997_loose';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'das1997_strict'} = 'data/alpha_polytopic_bitopic.das1997_strict';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'das2002'} = 'data/alpha_polytopic_bitopic.das2002';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'dastmfilter'} = 'data/alpha_polytopic_bitopic.dastmfilter';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'deltag'} = 'data/alpha_polytopic_bitopic.deltag';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'ensemble'} = 'data/alpha_polytopic_bitopic.ensemble';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'hmmtm'} = 'data/alpha_polytopic_bitopic.hmmtm';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'hmmtop2'} = 'data/alpha_polytopic_bitopic.hmmtop2';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'hmmtop_topcons_single'} = 'data/alpha_polytopic_bitopic.hmmtop_topcons_single';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'membrain'} = 'data/alpha_polytopic_bitopic.membrain';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'memsat3'} = 'data/alpha_polytopic_bitopic.memsat3';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'memsatsvm'} = 'data/alpha_polytopic_bitopic.memsatsvm';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'memsat_topcons_single'} = 'data/alpha_polytopic_bitopic.memsat_topcons_single';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'minnou'} = 'data/alpha_polytopic_bitopic.minnou';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'octopus'} = 'data/alpha_polytopic_bitopic.octopus';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'octopus_topcons'} = 'data/alpha_polytopic_bitopic.octopus_topcons';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'pepinfo_eisenberg_11_10'} = 'data/alpha_polytopic_bitopic.pepinfo_eisenberg_11_10';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'pepinfo_eisenberg_19_10'} = 'data/alpha_polytopic_bitopic.pepinfo_eisenberg_19_10';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'pepinfo_eisenberg_7_10'} = 'data/alpha_polytopic_bitopic.pepinfo_eisenberg_7_10';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'pepinfo_kyte_11_10'} = 'data/alpha_polytopic_bitopic.pepinfo_kyte_11_10';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'pepinfo_kyte_19_10'} = 'data/alpha_polytopic_bitopic.pepinfo_kyte_19_10';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'pepinfo_kyte_7_10'} = 'data/alpha_polytopic_bitopic.pepinfo_kyte_7_10';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'pepinfo_ohm_11_10'} = 'data/alpha_polytopic_bitopic.pepinfo_ohm_11_10';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'pepinfo_ohm_19_10'} = 'data/alpha_polytopic_bitopic.pepinfo_ohm_19_10';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'pepinfo_ohm_7_10'} = 'data/alpha_polytopic_bitopic.pepinfo_ohm_7_10';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'phdhtm_pbil'} = 'data/alpha_polytopic_bitopic.phdhtm_pbil';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'phdthtm_pbil'} = 'data/alpha_polytopic_bitopic.phdthtm_pbil';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'philius'} = 'data/alpha_polytopic_bitopic.philius';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'phobius'} = 'data/alpha_polytopic_bitopic.phobius';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'polyphobius'} = 'data/alpha_polytopic_bitopic.polyphobius';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'predtmr'} = 'data/alpha_polytopic_bitopic.predtmr';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'prodiv_topcons'} = 'data/alpha_polytopic_bitopic.prodiv_topcons';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'pro_topcons'} = 'data/alpha_polytopic_bitopic.pro_topcons';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'scampi'} = 'data/alpha_polytopic_bitopic.scampi';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'scampi_msa_topcons'} = 'data/alpha_polytopic_bitopic.scampi_msa_topcons';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'scampi_seq_topcons'} = 'data/alpha_polytopic_bitopic.scampi_seq_topcons';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'scampi_seq_topcons_single'} = 'data/alpha_polytopic_bitopic.scampi_seq_topcons_single';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'sosui'} = 'data/alpha_polytopic_bitopic.sosui';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'split4'} = 'data/alpha_polytopic_bitopic.split4';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'stmhmm_topcons_single'} = 'data/alpha_polytopic_bitopic.stmhmm_topcons_single';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'svmtm'} = 'data/alpha_polytopic_bitopic.svmtm';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'svmtop'} = 'data/alpha_polytopic_bitopic.svmtop';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'tmap'} = 'data/alpha_polytopic_bitopic.tmap';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'tmhmm2'} = 'data/alpha_polytopic_bitopic.tmhmm2';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'tmloop'} = 'data/alpha_polytopic_bitopic.tmloop';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'tmmod'} = 'data/alpha_polytopic_bitopic.tmmod';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'tmpred'} = 'data/alpha_polytopic_bitopic.tmpred';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'topcons'} = 'data/alpha_polytopic_bitopic.topcons';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'topcons_single'} = 'data/alpha_polytopic_bitopic.topcons_single';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'toppred2'} = 'data/alpha_polytopic_bitopic.toppred2';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'valpred2'} = 'data/alpha_polytopic_bitopic.valpred2';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'valpred'} = 'data/alpha_polytopic_bitopic.valpred';
	$data->{'file_name'}->{'alpha_polytopic_bitopic'}->{'wavetm'} = 'data/alpha_polytopic_bitopic.wavetm';

	$data->{'file_name'}->{'beta_barrel'}->{'das1997_loose'} = 'data/beta_barrel.das1997_loose';
	$data->{'file_name'}->{'beta_barrel'}->{'das1997_strict'} = 'data/beta_barrel.das1997_strict';
	$data->{'file_name'}->{'beta_barrel'}->{'das2002'} = 'data/beta_barrel.das2002';
	$data->{'file_name'}->{'beta_barrel'}->{'dastmfilter'} = 'data/beta_barrel.dastmfilter';
	$data->{'file_name'}->{'beta_barrel'}->{'deltag'} = 'data/beta_barrel.deltag';
	$data->{'file_name'}->{'beta_barrel'}->{'ensemble'} = 'data/beta_barrel.ensemble';
	$data->{'file_name'}->{'beta_barrel'}->{'hmmtm'} = 'data/beta_barrel.hmmtm';
	$data->{'file_name'}->{'beta_barrel'}->{'hmmtop2'} = 'data/beta_barrel.hmmtop2';
	$data->{'file_name'}->{'beta_barrel'}->{'hmmtop_topcons_single'} = 'data/beta_barrel.hmmtop_topcons_single';
	$data->{'file_name'}->{'beta_barrel'}->{'membrain'} = 'data/beta_barrel.membrain';
	$data->{'file_name'}->{'beta_barrel'}->{'memsat3'} = 'data/beta_barrel.memsat3';
	$data->{'file_name'}->{'beta_barrel'}->{'memsatsvm'} = 'data/beta_barrel.memsatsvm';
	$data->{'file_name'}->{'beta_barrel'}->{'memsat_topcons_single'} = 'data/beta_barrel.memsat_topcons_single';
	$data->{'file_name'}->{'beta_barrel'}->{'minnou'} = 'data/beta_barrel.minnou';
	$data->{'file_name'}->{'beta_barrel'}->{'octopus'} = 'data/beta_barrel.octopus';
	$data->{'file_name'}->{'beta_barrel'}->{'octopus_topcons'} = 'data/beta_barrel.octopus_topcons';
	$data->{'file_name'}->{'beta_barrel'}->{'pepinfo_eisenberg_11_10'} = 'data/beta_barrel.pepinfo_eisenberg_11_10';
	$data->{'file_name'}->{'beta_barrel'}->{'pepinfo_eisenberg_19_10'} = 'data/beta_barrel.pepinfo_eisenberg_19_10';
	$data->{'file_name'}->{'beta_barrel'}->{'pepinfo_eisenberg_7_10'} = 'data/beta_barrel.pepinfo_eisenberg_7_10';
	$data->{'file_name'}->{'beta_barrel'}->{'pepinfo_kyte_11_10'} = 'data/beta_barrel.pepinfo_kyte_11_10';
	$data->{'file_name'}->{'beta_barrel'}->{'pepinfo_kyte_19_10'} = 'data/beta_barrel.pepinfo_kyte_19_10';
	$data->{'file_name'}->{'beta_barrel'}->{'pepinfo_kyte_7_10'} = 'data/beta_barrel.pepinfo_kyte_7_10';
	$data->{'file_name'}->{'beta_barrel'}->{'pepinfo_ohm_11_10'} = 'data/beta_barrel.pepinfo_ohm_11_10';
	$data->{'file_name'}->{'beta_barrel'}->{'pepinfo_ohm_19_10'} = 'data/beta_barrel.pepinfo_ohm_19_10';
	$data->{'file_name'}->{'beta_barrel'}->{'pepinfo_ohm_7_10'} = 'data/beta_barrel.pepinfo_ohm_7_10';
	$data->{'file_name'}->{'beta_barrel'}->{'phdhtm_pbil'} = 'data/beta_barrel.phdhtm_pbil';
	$data->{'file_name'}->{'beta_barrel'}->{'phdthtm_pbil'} = 'data/beta_barrel.phdthtm_pbil';
	$data->{'file_name'}->{'beta_barrel'}->{'philius'} = 'data/beta_barrel.philius';
	$data->{'file_name'}->{'beta_barrel'}->{'phobius'} = 'data/beta_barrel.phobius';
	$data->{'file_name'}->{'beta_barrel'}->{'polyphobius'} = 'data/beta_barrel.polyphobius';
	$data->{'file_name'}->{'beta_barrel'}->{'predtmr'} = 'data/beta_barrel.predtmr';
	$data->{'file_name'}->{'beta_barrel'}->{'prodiv_topcons'} = 'data/beta_barrel.prodiv_topcons';
	$data->{'file_name'}->{'beta_barrel'}->{'pro_topcons'} = 'data/beta_barrel.pro_topcons';
	$data->{'file_name'}->{'beta_barrel'}->{'scampi'} = 'data/beta_barrel.scampi';
	$data->{'file_name'}->{'beta_barrel'}->{'scampi_msa_topcons'} = 'data/beta_barrel.scampi_msa_topcons';
	$data->{'file_name'}->{'beta_barrel'}->{'scampi_seq_topcons'} = 'data/beta_barrel.scampi_seq_topcons';
	$data->{'file_name'}->{'beta_barrel'}->{'scampi_seq_topcons_single'} = 'data/beta_barrel.scampi_seq_topcons_single';
	$data->{'file_name'}->{'beta_barrel'}->{'sosui'} = 'data/beta_barrel.sosui';
	$data->{'file_name'}->{'beta_barrel'}->{'split4'} = 'data/beta_barrel.split4';
	$data->{'file_name'}->{'beta_barrel'}->{'stmhmm_topcons_single'} = 'data/beta_barrel.stmhmm_topcons_single';
	$data->{'file_name'}->{'beta_barrel'}->{'svmtm'} = 'data/beta_barrel.svmtm';
	$data->{'file_name'}->{'beta_barrel'}->{'svmtop'} = 'data/beta_barrel.svmtop';
	$data->{'file_name'}->{'beta_barrel'}->{'tmap'} = 'data/beta_barrel.tmap';
	$data->{'file_name'}->{'beta_barrel'}->{'tmhmm2'} = 'data/beta_barrel.tmhmm2';
	$data->{'file_name'}->{'beta_barrel'}->{'tmloop'} = 'data/beta_barrel.tmloop';
	$data->{'file_name'}->{'beta_barrel'}->{'tmmod'} = 'data/beta_barrel.tmmod';
	$data->{'file_name'}->{'beta_barrel'}->{'tmpred'} = 'data/beta_barrel.tmpred';
	$data->{'file_name'}->{'beta_barrel'}->{'topcons'} = 'data/beta_barrel.topcons';
	$data->{'file_name'}->{'beta_barrel'}->{'topcons_single'} = 'data/beta_barrel.topcons_single';
	$data->{'file_name'}->{'beta_barrel'}->{'toppred2'} = 'data/beta_barrel.toppred2';
	$data->{'file_name'}->{'beta_barrel'}->{'valpred2'} = 'data/beta_barrel.valpred2';
	$data->{'file_name'}->{'beta_barrel'}->{'valpred'} = 'data/beta_barrel.valpred';
	$data->{'file_name'}->{'beta_barrel'}->{'wavetm'} = 'data/beta_barrel.wavetm';

	$data->{'file_name'}->{'soluble'}->{'das1997_loose'} = 'data/soluble.das1997_loose';
	$data->{'file_name'}->{'soluble'}->{'das1997_strict'} = 'data/soluble.das1997_strict';
	$data->{'file_name'}->{'soluble'}->{'das2002'} = 'data/soluble.das2002';
	$data->{'file_name'}->{'soluble'}->{'dastmfilter'} = 'data/soluble.dastmfilter';
	$data->{'file_name'}->{'soluble'}->{'deltag'} = 'data/soluble.deltag';
	$data->{'file_name'}->{'soluble'}->{'ensemble'} = 'data/soluble.ensemble';
	$data->{'file_name'}->{'soluble'}->{'hmmtm'} = 'data/soluble.hmmtm';
	$data->{'file_name'}->{'soluble'}->{'hmmtop2'} = 'data/soluble.hmmtop2';
	$data->{'file_name'}->{'soluble'}->{'hmmtop_topcons_single'} = 'data/soluble.hmmtop_topcons_single';
	$data->{'file_name'}->{'soluble'}->{'membrain'} = 'data/soluble.membrain';
	$data->{'file_name'}->{'soluble'}->{'memsat3'} = 'data/soluble.memsat3';
	$data->{'file_name'}->{'soluble'}->{'memsatsvm'} = 'data/soluble.memsatsvm';
	$data->{'file_name'}->{'soluble'}->{'memsat_topcons_single'} = 'data/soluble.memsat_topcons_single';
	$data->{'file_name'}->{'soluble'}->{'minnou'} = 'data/soluble.minnou';
	$data->{'file_name'}->{'soluble'}->{'octopus'} = 'data/soluble.octopus';
	$data->{'file_name'}->{'soluble'}->{'octopus_topcons'} = 'data/soluble.octopus_topcons';
	$data->{'file_name'}->{'soluble'}->{'pepinfo_eisenberg_11_10'} = 'data/soluble.pepinfo_eisenberg_11_10';
	$data->{'file_name'}->{'soluble'}->{'pepinfo_eisenberg_19_10'} = 'data/soluble.pepinfo_eisenberg_19_10';
	$data->{'file_name'}->{'soluble'}->{'pepinfo_eisenberg_7_10'} = 'data/soluble.pepinfo_eisenberg_7_10';
	$data->{'file_name'}->{'soluble'}->{'pepinfo_kyte_11_10'} = 'data/soluble.pepinfo_kyte_11_10';
	$data->{'file_name'}->{'soluble'}->{'pepinfo_kyte_19_10'} = 'data/soluble.pepinfo_kyte_19_10';
	$data->{'file_name'}->{'soluble'}->{'pepinfo_kyte_7_10'} = 'data/soluble.pepinfo_kyte_7_10';
	$data->{'file_name'}->{'soluble'}->{'pepinfo_ohm_11_10'} = 'data/soluble.pepinfo_ohm_11_10';
	$data->{'file_name'}->{'soluble'}->{'pepinfo_ohm_19_10'} = 'data/soluble.pepinfo_ohm_19_10';
	$data->{'file_name'}->{'soluble'}->{'pepinfo_ohm_7_10'} = 'data/soluble.pepinfo_ohm_7_10';
	$data->{'file_name'}->{'soluble'}->{'phdhtm_pbil'} = 'data/soluble.phdhtm_pbil';
	$data->{'file_name'}->{'soluble'}->{'phdthtm_pbil'} = 'data/soluble.phdthtm_pbil';
	$data->{'file_name'}->{'soluble'}->{'philius'} = 'data/soluble.philius';
	$data->{'file_name'}->{'soluble'}->{'phobius'} = 'data/soluble.phobius';
	$data->{'file_name'}->{'soluble'}->{'polyphobius'} = 'data/soluble.polyphobius';
	$data->{'file_name'}->{'soluble'}->{'predtmr'} = 'data/soluble.predtmr';
	$data->{'file_name'}->{'soluble'}->{'prodiv_topcons'} = 'data/soluble.prodiv_topcons';
	$data->{'file_name'}->{'soluble'}->{'pro_topcons'} = 'data/soluble.pro_topcons';
	$data->{'file_name'}->{'soluble'}->{'scampi'} = 'data/soluble.scampi';
	$data->{'file_name'}->{'soluble'}->{'scampi_msa_topcons'} = 'data/soluble.scampi_msa_topcons';
	$data->{'file_name'}->{'soluble'}->{'scampi_seq_topcons'} = 'data/soluble.scampi_seq_topcons';
	$data->{'file_name'}->{'soluble'}->{'scampi_seq_topcons_single'} = 'data/soluble.scampi_seq_topcons_single';
	$data->{'file_name'}->{'soluble'}->{'sosui'} = 'data/soluble.sosui';
	$data->{'file_name'}->{'soluble'}->{'split4'} = 'data/soluble.split4';
	$data->{'file_name'}->{'soluble'}->{'stmhmm_topcons_single'} = 'data/soluble.stmhmm_topcons_single';
	$data->{'file_name'}->{'soluble'}->{'svmtm'} = 'data/soluble.svmtm';
	$data->{'file_name'}->{'soluble'}->{'svmtop'} = 'data/soluble.svmtop';
	$data->{'file_name'}->{'soluble'}->{'tmap'} = 'data/soluble.tmap';
	$data->{'file_name'}->{'soluble'}->{'tmhmm2'} = 'data/soluble.tmhmm2';
	$data->{'file_name'}->{'soluble'}->{'tmloop'} = 'data/soluble.tmloop';
	$data->{'file_name'}->{'soluble'}->{'tmmod'} = 'data/soluble.tmmod';
	$data->{'file_name'}->{'soluble'}->{'tmpred'} = 'data/soluble.tmpred';
	$data->{'file_name'}->{'soluble'}->{'topcons'} = 'data/soluble.topcons';
	$data->{'file_name'}->{'soluble'}->{'topcons_single'} = 'data/soluble.topcons_single';
	$data->{'file_name'}->{'soluble'}->{'toppred2'} = 'data/soluble.toppred2';
	$data->{'file_name'}->{'soluble'}->{'valpred2'} = 'data/soluble.valpred2';
	$data->{'file_name'}->{'soluble'}->{'valpred'} = 'data/soluble.valpred';
	$data->{'file_name'}->{'soluble'}->{'wavetm'} = 'data/soluble.wavetm';

	#----------

	my %data_method_short_name;
	$data->{'method_short_name'} = \%data_method_short_name;
	$data->{'method_short_name'}->{'dssp'} = 'DSSP';
	$data->{'method_short_name'}->{'opm_membrane_near'} = 'OPM-Mnear';
	$data->{'method_short_name'}->{'opm_posn'} = 'OPMtoplgy';
	$data->{'method_short_name'}->{'opm_seq'} = 'PDBfile';
	$data->{'method_short_name'}->{'opm_adj_tm_segments'} = 'OPM-adj';
	$data->{'method_short_name'}->{'opm_adj_tm_segments_topology'} = 'OPMadjTOP';
	$data->{'method_short_name'}->{'opm_adj_tmh_tmbb_segments'} = 'OPM-D-adj';
	$data->{'method_short_name'}->{'opm_tm_segments'} = 'OPM';
	$data->{'method_short_name'}->{'opm_tm_segments_topology'} = 'OPM___TOP';
	$data->{'method_short_name'}->{'opm_tmh_tmbb_segments'} = 'OPM-DSSP';
	$data->{'method_short_name'}->{'opm_tm_subunits'} = 'OPMsubunt';
	$data->{'method_short_name'}->{'pdb_seq'} = 'SEQUENCE';
	$data->{'method_short_name'}->{'pdbtm_segments'} = 'PDBTM';
	$data->{'method_short_name'}->{'pdbtm_segments_topology'} = 'PDBTM_TOP';
	$data->{'method_short_name'}->{'stride'} = 'STRIDE';
	$data->{'method_short_name'}->{'stride_seq'} = 'PDBfile2';
	#----------
	$data->{'method_short_name'}->{'you'} = 'YOU';
	#----------
	$data->{'method_short_name'}->{'dastmfilter'} = 'DASTMfilt';
	$data->{'method_short_name'}->{'das2002'} = 'DAS2002';
	$data->{'method_short_name'}->{'das1997_loose'} = 'DAS1997l';
	$data->{'method_short_name'}->{'das1997_strict'} = 'DAS1997s';
	$data->{'method_short_name'}->{'deltag'} = 'deltaG';
	$data->{'method_short_name'}->{'ensemble'} = 'ENSEMBLE';
	$data->{'method_short_name'}->{'hmmtm'} = 'HMM-TM';
	$data->{'method_short_name'}->{'hmmtop2'} = 'HMMTOP2';
	$data->{'method_short_name'}->{'hmmtop_topcons_single'} = 'hmmtopS';
	$data->{'method_short_name'}->{'membrain'} = 'MemBrain';
	$data->{'method_short_name'}->{'memsatsvm'} = 'MEMSATSVM';
	$data->{'method_short_name'}->{'memsat3'} = 'MEMSAT3';
	$data->{'method_short_name'}->{'memsat_topcons_single'} = 'memsatS';
	$data->{'method_short_name'}->{'minnou'} = 'MINNOU';
	$data->{'method_short_name'}->{'octopus'} = 'OCTOPUS';
	$data->{'method_short_name'}->{'octopus_topcons'} = 'octopusT';
	$data->{'method_short_name'}->{'phdhtm_pbil'} = 'PHDhtm';
	$data->{'method_short_name'}->{'phdthtm_pbil'} = 'PHDThtm';
	$data->{'method_short_name'}->{'philius'} = 'Philius';
	$data->{'method_short_name'}->{'phobius'} = 'Phobius';
	$data->{'method_short_name'}->{'polyphobius'} = 'PolyPhobs';
	$data->{'method_short_name'}->{'predtmr'} = 'PRED-TMR';
	$data->{'method_short_name'}->{'pro_topcons'} = 'proT';
	$data->{'method_short_name'}->{'prodiv_topcons'} = 'prodivT';
	$data->{'method_short_name'}->{'scampi'} = 'SCAMPI';
	$data->{'method_short_name'}->{'scampi_seq_topcons'} = 'SCAMPIsqT';
	$data->{'method_short_name'}->{'scampi_msa_topcons'} = 'SCAMPImaT';
	$data->{'method_short_name'}->{'scampi_seq_topcons_single'} = 'SCAMPIsqS';
	$data->{'method_short_name'}->{'sosui'} = 'SOSUI';
	$data->{'method_short_name'}->{'split4'} = 'SPLIT4';
	$data->{'method_short_name'}->{'svmtm'} = 'SVMtm';
	$data->{'method_short_name'}->{'svmtop'} = 'SVMtop';
	$data->{'method_short_name'}->{'tmap'} = 'TMAP';
	$data->{'method_short_name'}->{'tmhmm2'} = 'TMHMM2';
	$data->{'method_short_name'}->{'stmhmm_topcons_single'} = 'stmhmmS';
	$data->{'method_short_name'}->{'tmloop'} = 'TMLOOP';
	$data->{'method_short_name'}->{'tmmod'} = 'TMMOD';
	$data->{'method_short_name'}->{'tmpred'} = 'TMpred';
	$data->{'method_short_name'}->{'topcons'} = 'TOPCONS';
	$data->{'method_short_name'}->{'topcons_single'} = 'TOPCONSs';
	$data->{'method_short_name'}->{'toppred2'} = 'TOPPRED2';
	$data->{'method_short_name'}->{'valpred'} = 'VALPRED';
	$data->{'method_short_name'}->{'valpred2'} = 'VALPRED2';
	$data->{'method_short_name'}->{'wavetm'} = 'waveTM';
	$data->{'method_short_name'}->{'pepinfo_eisenberg_7_10'} = 'Eisen(7)';
	$data->{'method_short_name'}->{'pepinfo_eisenberg_11_10'} = 'Eisen(11)';
	$data->{'method_short_name'}->{'pepinfo_eisenberg_19_10'} = 'Eisen(19)';
	$data->{'method_short_name'}->{'pepinfo_kyte_7_10'} = 'KyteD(7)';
	$data->{'method_short_name'}->{'pepinfo_kyte_11_10'} = 'KyteD(11)';
	$data->{'method_short_name'}->{'pepinfo_kyte_19_10'} = 'KyteD(19)';
	$data->{'method_short_name'}->{'pepinfo_ohm_7_10'} = 'OHM(7)';
	$data->{'method_short_name'}->{'pepinfo_ohm_11_10'} = 'OHM(11)';
	$data->{'method_short_name'}->{'pepinfo_ohm_19_10'} = 'OHM(19)';
	$data->{'method_short_name'}->{'you'} = 'YOU';
	#----------
	$data->{'method_long_name'}->{'dastmfilter'} = 'DAS-TMfilter';
	$data->{'method_long_name'}->{'das2002'} = 'DAS2002';
	$data->{'method_long_name'}->{'das1997_loose'} = 'DAS1997 (loose)';
	$data->{'method_long_name'}->{'das1997_strict'} = 'DAS1997 (strict)';
	$data->{'method_long_name'}->{'deltag'} = 'deltaG';
	$data->{'method_long_name'}->{'ensemble'} = 'ENSEMBLE (in MemPype)';
	$data->{'method_long_name'}->{'hmmtm'} = 'HMM-TM';
	$data->{'method_long_name'}->{'hmmtop2'} = 'HMMTOP2';
	$data->{'method_long_name'}->{'hmmtop_topcons_single'} = 'HMMTOP (in TOPCONS-single)';
	$data->{'method_long_name'}->{'membrain'} = 'MemBrain';
	$data->{'method_long_name'}->{'memsatsvm'} = 'MEMSAT-SVM';
	$data->{'method_long_name'}->{'memsat3'} = 'MEMSAT3';
	$data->{'method_long_name'}->{'memsat_topcons_single'} = 'MEMSAT (in TOPCONS-single)';
	$data->{'method_long_name'}->{'minnou'} = 'MINNOU';
	$data->{'method_long_name'}->{'octopus'} = 'OCTOPUS';
	$data->{'method_long_name'}->{'octopus_topcons'} = 'OCTOPUS (in TOPCONS)';
	$data->{'method_long_name'}->{'phdhtm_pbil'} = 'PHDhtm (at PBIL)';
	$data->{'method_long_name'}->{'phdthtm_pbil'} = 'PHDThtm (at PBIL)';
	$data->{'method_long_name'}->{'philius'} = 'Philius';
	$data->{'method_long_name'}->{'phobius'} = 'Phobius';
	$data->{'method_long_name'}->{'polyphobius'} = 'PolyPhobius';
	$data->{'method_long_name'}->{'predtmr'} = 'PRED-TMR';
	$data->{'method_long_name'}->{'pro_topcons'} = 'PRO-TMHMM (in TOPCONS)';
	$data->{'method_long_name'}->{'prodiv_topcons'} = 'PRODIV-TMHMM (in TOPCONS)';
	$data->{'method_long_name'}->{'scampi'} = 'SCAMPI';
	$data->{'method_long_name'}->{'scampi_seq_topcons'} = 'SCAMPI-sequence (in TOPCONS)';
	$data->{'method_long_name'}->{'scampi_msa_topcons'} = 'SCAMPI-multi (in TOPCONS)';
	$data->{'method_long_name'}->{'scampi_seq_topcons_single'} = 'SCAMPI-sequence (in TOPCONS-single)';
	$data->{'method_long_name'}->{'sosui'} = 'SOSUI';
	$data->{'method_long_name'}->{'split4'} = 'SPLIT4';
	$data->{'method_long_name'}->{'svmtm'} = 'SVMtm';
	$data->{'method_long_name'}->{'svmtop'} = 'SVMtop';
	$data->{'method_long_name'}->{'tmap'} = 'TMAP';
	$data->{'method_long_name'}->{'tmhmm2'} = 'TMHMM2';
	$data->{'method_long_name'}->{'stmhmm_topcons_single'} = 'S-TMHMM (in TOPCONS-single)';
	$data->{'method_long_name'}->{'tmloop'} = 'TMLOOP';
	$data->{'method_long_name'}->{'tmmod'} = 'TMMOD';
	$data->{'method_long_name'}->{'tmpred'} = 'TMpred';
	$data->{'method_long_name'}->{'topcons'} = 'TOPCONS';
	$data->{'method_long_name'}->{'topcons_single'} = 'TOPCONS-single';
	$data->{'method_long_name'}->{'toppred2'} = 'TOPPRED2';
	$data->{'method_long_name'}->{'valpred'} = 'VALPRED';
	$data->{'method_long_name'}->{'valpred2'} = 'VALPRED2';
	$data->{'method_long_name'}->{'wavetm'} = 'waveTM';
	$data->{'method_long_name'}->{'pepinfo_eisenberg_7_10'} = 'Eisen(7,10)';
	$data->{'method_long_name'}->{'pepinfo_eisenberg_11_10'} = 'Eisen(11,10)';
	$data->{'method_long_name'}->{'pepinfo_eisenberg_19_10'} = 'Eisen(19,10)';
	$data->{'method_long_name'}->{'pepinfo_kyte_7_10'} = 'KyteD(7,10)';
	$data->{'method_long_name'}->{'pepinfo_kyte_11_10'} = 'KyteD(11,10)';
	$data->{'method_long_name'}->{'pepinfo_kyte_19_10'} = 'KyteD(19,10)';
	$data->{'method_long_name'}->{'pepinfo_ohm_7_10'} = 'OHM(7,10)';
	$data->{'method_long_name'}->{'pepinfo_ohm_11_10'} = 'OHM(11,10)';
	$data->{'method_long_name'}->{'pepinfo_ohm_19_10'} = 'OHM(19,10)';
	$data->{'method_long_name'}->{'you'} = 'YOU';
	#----------

	my @output_prediction_method;
	$output->{'prediction_method'} = \@output_prediction_method;
	my %output_prediction_sequence;
	$output->{'prediction_sequence'} = \%output_prediction_sequence;

	$output->{'num_reference_seqs'} = 0;
	$output->{'num_reference_membrane_helices'} = 0;
	$output->{'num_user_seqs'} = 0;
	$output->{'num_user_membrane_helices'} = 0;

	$global->{'benchmark_standard'} = '';
	$global->{'show_intermediate'} = 0;
	$global->{'show_topology_files'} = 0;
	$output->{'display_method'}->{'opm_adj_tm_segments'} = 0;
	$output->{'display_method'}->{'opm_adj_tm_segments_topology'} = 0;
	$output->{'display_method'}->{'opm_tm_segments'} = 0;
	$output->{'display_method'}->{'opm_tm_segments_topology'} = 0;
	$output->{'display_method'}->{'pdbtm_segments'} = 0;
	$output->{'display_method'}->{'pdbtm_segments_topology'} = 0;
}

#===================================================================================================
#	make sure the user-entered a positive number (can be decimals, can be zero)
#===================================================================================================

sub is_valid_positive_number_from_zero {

	my $input_number = shift;
	my $return = 1;

	if ($input_number !~ /^\d+\.?\d*$/) {
		$return = 0;
	} else {
		if ($input_number < 0) {
			$return = 0;
		}
	}

	return $return;
}

#===================================================================================================
#	make sure the user-entered a positive or negative number (can be decimals, can be zero)
#===================================================================================================

sub is_valid_number {

	my $input_number = shift;
	my $return = 1;

	if ($input_number !~ /^-?\d+\.?\d*$/) {
		$return = 0;
	}

	return $return;
}

#===================================================================================================
#	make sure the user-entered number is 0, 1, 2, ...
#===================================================================================================

sub is_valid_positive_integer_from_zero {

	my $input_number = shift;
	my $return = 1;

	if ($input_number !~ /^\d+$/) {
		$return = 0;
	} else {
		if ($input_number < 0) {
			$return = 0;
		}
	}

	return $return;
}

#===================================================================================================
#	make sure the user-entered number is 1, 2, 3, ...
#===================================================================================================

sub is_valid_positive_integer_above_zero {

	my $input_number = shift;
	my $return = 1;

	if ($input_number !~ /^\d+$/) {
		$return = 0;
	} else {
		if ($input_number < 1) {
			$return = 0;
		}
	}

	return $return;
}

#===================================================================================================
#	Perl trim function to remove whitespace from the start and end of the string
#===================================================================================================

sub trim {
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}
