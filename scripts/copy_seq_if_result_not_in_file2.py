#!/usr/bin/env python

# python copy_seq_if_result_not_in_file2.py alpha_polytopic_bitopic_noblanks.fasta alpha_polytopic_bitopic_noblanks.fasta.memsatsvm_link_1 alpha_polytopic_bitopic_noblanks.fasta.memsatsvm2

# python copy_seq_if_result_not_in_file2.py [input-file-1] [input-file-2] [output-file]

import sys
import re


# global variables
save_id = []


######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	output_file_name = sys.argv[3]

	infile2 = open( input_file_2, "r" )
	inlines2 = infile2.readlines()
	infile2.close()

	for inline in inlines2:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_result = bits[1]
			if (this_result != '-'):
				save_id.append( this_id )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	# parse the sequences because they are over one line

	outfile = open( output_file_name, "w" )

	output_this_id = 0
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):

			if (inline[0:1] == '>'):
				curr_id = inline[1:7]

				this_id_is_in_file2 = 0
				for this_save_id in save_id:
					if (this_save_id == curr_id):
						this_id_is_in_file2 = 1

				if (this_id_is_in_file2 == 1):
					output_this_id = 0
				else:
					output_this_id = 1

			if (output_this_id == 1):
				output_line = inline + "\r\n"
				outfile.write( output_line )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

