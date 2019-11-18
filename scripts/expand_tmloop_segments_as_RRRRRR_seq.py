#!/usr/bin/env python

# python expand_tmloop_segments_as_RRRRRR_seq.py alpha_polytopic_biotopic_noblanks.tmloop.adjusted_for_removed_blanks alpha_polytopic_biotopic.fasta alpha_polytopic_biotopic.tmloop

# python expand_tmloop_segments_as_RRRRRR_seq.py [input-file-1] [input-file-2] [output-file]

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
	output_file_name = sys.argv[3]

	infile2 = open( input_file_2, "r" )
	inlines2 = infile2.readlines()
	infile2.close()

	curr_id = ''
	curr_seq = ''
	for inline in inlines2:
		inline = inline.strip()
		if (inline != ''):
			if (inline[0:1] == '>'):
				if (curr_id != ''):
					save_id.append( curr_id )
					save_seq_length.append( len(curr_seq) )
				curr_id = inline[1:7]
				curr_seq = ''
			else:
				curr_seq = curr_seq + inline
	if (curr_id != ''):
		save_id.append( curr_id )
		save_seq_length.append( len(curr_seq) )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	outfile = open( output_file_name, "w" )

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_segments = bits[1]
			this_seq_length = 0

			i = 0
			for this_save_id in save_id:
				if (this_save_id == this_id):
					this_seq_length = save_seq_length[i]
				i = i + 1
			if (this_seq_length == 0):
				print this_id, ': ERROR : Sequence length is zero'

			new_seq_array = []
			for i in range( 0, this_seq_length ):
				new_seq_array.append( '_' )
			
			if (this_segments != '-'):
				bits = this_segments.rsplit(',')
				for bit in bits:
					bit = bit.strip()
					if (bit != ''):
						bits2 = bit.rsplit('-')
						this_from = int(bits2[0])
						this_to = int(bits2[1])
						for i in range( (this_from - 1), this_to ):
							if (i < this_seq_length):
								new_seq_array[i] = 'R'

			new_seq = ''
			for i in range( 0, this_seq_length ):
				new_seq = new_seq + new_seq_array[i]
			if ((new_seq == '') or (new_seq == '_')):
				new_seq = '-'

			output_line = this_id + ":::::" + new_seq + ":::::" + "\r\n"
			outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

