#!/usr/bin/env python

# python parse_PHDhtm_PBIL_results.py alpha_polytopic_bitopic_noblanks.fasta.PHDhtm_PBIL_output alpha_polytopic_bitopic_noblanks.fasta alpha_polytopic_bitopic_noblanks.PHDhtm_PBIL_results

# python parse_PHDhtm_PBIL_results.py [input-file] [input-file-2] [output-file-hmmtm]

import sys
import re


# global variables
save_id = []
save_seq_length = []


######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	output_file = sys.argv[3]

	# read in the lengths of each sequence, contained in input file 2

	infile2 = open( input_file_2, "r" )
	inlines2 = infile2.readlines()
	infile2.close()

	in_a_seq = 0
	curr_id = ''
	curr_seq = ''
	for inline in inlines2:
		inline = inline.strip()
		if (inline != ''):
			if (inline[0:1] == '>'):
				if (in_a_seq == 1):
					save_id.append( curr_id )
					save_seq_length.append( len(curr_seq) )
				in_a_seq = 1
				curr_id = inline[1:7]
				curr_seq = ''
			else:
				curr_seq = curr_seq + inline
	if (in_a_seq == 1):
		save_id.append( curr_id )
		save_seq_length.append( len(curr_seq) )

	# process input file 1 that contains the results to parse

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	outfile = open( output_file, "w" )

	curr_id = ''
	curr_string_array = []
	curr_seq_length = 0
	helix_count = 0
	curr_string = ''
	i = 0
	for inline in inlines1:

		is_id_line = 0
		if ( len(inline) >= 11):
			if (inline[6:11] == ':::::'):
				is_id_line = 1

		# 2xzb_A:::::
		if (is_id_line == 1):

			# finish previous id
			#if (curr_id != ''):
			#	for j in range( 0, curr_seq_length ):
			#		curr_string = curr_string + curr_string_array[j]
			#	output_line = curr_id + ':::::' + curr_string + ':::::' + "\r\n"
			#	outfile.write( output_line )

			# initialise for new id
			curr_id = inline[0:6]
			curr_string_array = []
			curr_seq_length = 0
			curr_string = ''
			found_id = 0
			j = 0
			for this_save_id in save_id:
				if (this_save_id == curr_id):
					found_id = 1
					curr_seq_length = save_seq_length[j]
				j = j + 1
			if (found_id == 1):
				for j in range( 0, curr_seq_length ):
					curr_string_array.append( '_' )
			else:
					print curr_id, 'did not find its sequence length'

		else:

			if ( inline.find('Protein Sci. 1995 Mar;4(3):521-33.</FONT><PRE><CODE>') != -1 ):

				#					Protein Sci. 1995 Mar;4(3):521-33.</FONT><PRE><CODE>
				#        10        20        30        40        50        60        70        80        90       100       110       120       130       140       150
				#         |         |         |         |         |         |         |         |         |         |         |         |         |         |         |
				#DIAGLTPCKESKQFAKREKQALKKLQASLKLYADDSAPALAIKATMEKTKKRFDNYGKYGLLCGSDGLPHLIVSGDQRHWGEFITPGILFLYIAGWIGWVGRSYLIAIRDEKKPTQKEIIIDVPLASSLLFRGFSWPVAAYRELLNGELVDNNF
				#                                                                                    <FONT COLOR=red>HHHHHHHHHHHHHHHH</FONT>                         <FONT COLOR=red>HHHHHHHHHHHH</FONT>                 

				this_row = inlines1[ i + 4 ]
				this_row = this_row.replace( '<FONT COLOR=red>', '' )
				this_row = this_row.replace( '</FONT>', '' )
				this_row = this_row.replace( ' ', '_' )
				this_row = this_row.replace( "\r", '' )
				this_row = this_row.replace( "\n", '' )
				curr_string = this_row
				output_line = curr_id + ':::::' + curr_string + ':::::' + "\r\n"
				outfile.write( output_line )

		i = i + 1

	# finish last id
	#if (curr_id != ''):
	#	for j in range( 0, curr_seq_length ):
	#		curr_string = curr_string + curr_string_array[j]
	#	output_line = curr_id + ':::::' + curr_string + ':::::' + "\r\n"
	#	outfile.write( output_line )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()


# input_file :

#2o01_F:::::
#					Protein Sci. 1995 Mar;4(3):521-33.</FONT><PRE><CODE>
#        10        20        30        40        50        60        70        80        90       100       110       120       130       140       150
#         |         |         |         |         |         |         |         |         |         |         |         |         |         |         |
#DIAGLTPCKESKQFAKREKQALKKLQASLKLYADDSAPALAIKATMEKTKKRFDNYGKYGLLCGSDGLPHLIVSGDQRHWGEFITPGILFLYIAGWIGWVGRSYLIAIRDEKKPTQKEIIIDVPLASSLLFRGFSWPVAAYRELLNGELVDNNF
#                                                                                    <FONT COLOR=red>HHHHHHHHHHHHHHHH</FONT>                         <FONT COLOR=red>HHHHHHHHHHHH</FONT>                 
#
#
#
#</CODE></PRE>
#<P>
#Prediction result file (text): [<A HREF=/tmp/7317801ceb27.phdhtm>PHD</A>]
#<BR>
#Intermediate result file (text): [<A HREF=/tmp/7317801ceb27.blastpsecpred>BLASTP on NRPROT</A>] [<A HREF=/tmp/7317801ceb27.cluoutphdpred>CLUSTALW</A>]
#</P>
#3rfr_F:::::

