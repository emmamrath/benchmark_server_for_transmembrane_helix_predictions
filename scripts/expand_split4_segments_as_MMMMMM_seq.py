#!/usr/bin/env python

# python expand_split4_segments_as_MMMMMM_seq.py alpha_polytopic_bitopic_noblanks.fasta.split4.adjusted_for_removed_blanks alpha_polytopic_bitopic.split4

# python expand_split4_segments_as_MMMMMM_seq.py [input-file] [output-file]

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
			this_topology = bits[3]
			output_seq = '-'
			output_seq_array = []
			prev_from = -1
			prev_to = -1
			topology_char = this_topology
			if (this_segments != '-'):
				output_seq = ''
				for i in range( 0, len(this_seq) ):
					output_seq_array.append( '_' )
				found_start_residue = 0
				for i in range( 0, len(this_seq) ):
					if (found_start_residue == 0):
						j = i + 1
						this_char = this_seq[i:j]
						if (this_char != '_'):
							prev_to = i + 1 - 1
							found_start_residue = 1
				bits2 = this_segments.rsplit( ',' )
				for this_segment in bits2:
					if (this_segment != ''):
						this_segment = this_segment.replace( '--', '-' )
						if (this_segment[0:1] == '-'):
							this_segment = this_segment[1:]
						bits3 = this_segment.rsplit( '-' )
						this_from = bits3[0]
						this_to = bits3[1]
						this_from = int(this_from)
						this_to = int(this_to)
						for j in range( prev_to, (this_from - 1) ):
							output_seq_array[j] = topology_char
						for j in range( (this_from - 1), this_to ):
							output_seq_array[j] = 'M'
						prev_from = this_from
						prev_to = this_to
						if (topology_char == 'i'):
							topology_char = 'o'
						elif (topology_char == 'o'):
							topology_char = 'i'
				for i in range( prev_to, len(this_seq) ):
					output_seq_array[i] = topology_char
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

