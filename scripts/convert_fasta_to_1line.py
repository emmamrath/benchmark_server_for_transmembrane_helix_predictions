#!/usr/bin/env python

# python convert_fasta_to_1line.py alpha_polytopic_bitopic_noblanks.fasta alpha_polytopic_bitopic_noblanks.1line

# python convert_fasta_to_1line.py [input-file] [output-file]

import sys
import re


# global variables
save_id = []
save_seq = []

######################################################
def is_valid_resn( this_char ):
	is_valid = 0
	for this_resn in resn_list_1char:
		if (this_resn == this_char):
			is_valid = 1
	return is_valid

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	output_file_name = sys.argv[2]

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	# parse the sequences because they are over one line

	in_a_seq = 0
	curr_id = ''
	curr_seq = ''
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			if (inline[0:1] == '>'):
				if (in_a_seq == 1):
					save_id.append( curr_id )
					save_seq.append( curr_seq )
				in_a_seq = 1
				curr_id = inline[1:7]
				curr_seq = ''
			else:
				curr_seq = curr_seq + inline
	if (in_a_seq == 1):
		save_id.append( curr_id )
		save_seq.append( curr_seq )

	outfile = open( output_file_name, "w" )
	i = 0
	for this_id in save_id:
		this_seq = save_seq[i]
		this_fastaid = '>' + this_id
		output_line = this_id + ':::::' + this_fastaid + ':::::' + this_seq + ':::::' + "\r\n"
		outfile.write( output_line )
		i = i + 1
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

