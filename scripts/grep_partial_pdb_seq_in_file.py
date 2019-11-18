#!/usr/bin/env python

# python ../grep_partial_pdb_seq_in_file.py list11_not_in_mpstruc_pass1_nor_pass2.pdb_seq one2_mpstruc_category.id.pdb_seq list11_not_in_mpstruc_pass1_nor_pass2.possible_ids

# python grep_partial_pdb_seq_in_file.py [input-file-1] [input-file-2] [output-file]

import sys
import re
import os
import subprocess

# global variables
save_id = []
save_seq = []


######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	output_file_name = sys.argv[3]

	#infile2 = open( input_file_2, "r" )
	#inlines2 = infile2.readlines()
	#infile2.close()

	#for inline in inlines2:
	#	inline = inline.strip()
	#	if (inline != ''):
	#		bits = inline.rsplit( ':::::' )
	#		this_id = bits[0]
	#		this_seq = bits[1]

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	outfile = open( output_file_name, "w" )

	global search_seqs
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			search_seq = this_seq
			this_len = len(search_seq)
			one_tenth = int(this_len / 10)
			one_fifth = int(this_len / 5)
			#search_seq = search_seq[ one_tenth : (this_len - one_tenth) ]
			search_seq = search_seq[ one_fifth : (this_len - one_fifth) ]
			#if ( this_len > 50 ):
			#	search_seq = search_seq[ 10 : (this_len - 10) ]

			#this_result = subprocess.call( [ "grep", search_seq, input_file_2 ] )			

			process = subprocess.Popen( [ "grep", search_seq, input_file_2 ], shell=False, stdout=subprocess.PIPE )
			this_result = process.communicate()

			if (this_result[0] != ''):
				grep_lines = this_result[0]
				grep_results = grep_lines.rsplit( "\r\n" )
				#for this_grep_result in grep_results:
				first_grep_result = grep_results[0]
				bits2 = first_grep_result.rsplit( ':::::' )
				grep_id = bits2[0]

				output_line = this_id + ':::::' + grep_id + ':::::' + "\r\n"
				outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

