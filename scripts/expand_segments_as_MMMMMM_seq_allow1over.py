#!/usr/bin/env python

# python expand_segments_as_MMMMMM_seq_allow1over.py alpha_polytopic_bitopic_noblanks.1line.das2002.adjusted_for_removed_blanks alpha_polytopic_bitopic.das2002

# python expand_segments_as_MMMMMM_seq_allow1over.py [input-file] [output-file]

import sys
import re


# global variables


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

	outfile = open( output_file_name, "w" )

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			this_segments = bits[2]
			output_seq = '-'
			output_seq_array = []
			if (this_segments != '-'):
				output_seq = ''
				for i in range( 0, len(this_seq) ):
					output_seq_array.append( '_' )
				length_seq = len(this_seq)
				length_seq_plus_1 = len(this_seq) + 1
				bits2 = this_segments.rsplit( ',' )
				for this_segment in bits2:
					if (this_segment != ''):
						bits3 = this_segment.rsplit( '-' )
						this_from = bits3[0]
						this_to = bits3[1]
						this_from = int(this_from)
						this_to = int(this_to)
						if (this_to == length_seq_plus_1):
							this_to = length_seq
						for j in range( (this_from - 1), this_to ):
							output_seq_array[j] = 'M'
				for i in range( 0, len(this_seq) ):
					output_seq = output_seq + output_seq_array[i]
			if ((output_seq == '_') or (output_seq == '')):
				output_seq = '-'
			output_line = this_id + ":::::" + output_seq + ":::::" + "\r\n"
			outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

