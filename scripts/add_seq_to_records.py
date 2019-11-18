#!/usr/bin/env python

# python add_seq_to_records.py alpha_polytopic_bitopic_noblanks.fasta.membrain.adjusted_for_removed_blanks alpha_polytopic_bitopic.pdb_seq alpha_polytopic_bitopic_noblanks.fasta.membrain.adjusted_for_removed_blanks.with_seq

# python add_seq_to_records.py [input-file] [input-file-2] [output-file]

import sys
import re


# global variables
save_id = []
save_seq = []


######################################################
def find_seq( this_id ):
	this_seq = ''
	i = 0
	for this_save_id in save_id:
		if (this_save_id == this_id):
			this_seq = save_seq[i]
		i = i + 1
	return this_seq

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
			this_seq = bits[1]
			save_id.append( this_id )
			save_seq.append( this_seq )

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
			this_seq = find_seq( this_id )
			if (this_seq == ''):
				print this_id, ': could not find a record for the sequence.'
			else:
				output_line = this_id + ':::::' + this_seq + ':::::' + this_segments + ':::::' + "\r\n"
				outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

