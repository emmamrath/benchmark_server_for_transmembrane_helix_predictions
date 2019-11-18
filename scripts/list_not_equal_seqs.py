#!/usr/bin/env python

# python ../list_not_equal_seqs.py alpha_polytopic_bitopic_peptide.opm_seq alpha_polytopic_bitopic_peptide.pdb_seq alpha_polytopic_bitopic_peptide.opm_seq_not_equal_pdb_seq.id

# python list_not_equal_seqs.py [input-file-1] [input-file-2] [output-file]

import sys
import re


# global variables
save_id = []
save_seq = []


######################################################
def return_seq_for_id( this_id ):

	found_seq = ''
	i = 0
	for this_save_id in save_id:
		if (this_save_id == this_id):
			found_seq = save_seq[i]
		i = i + 1
	return found_seq

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	output_file_name = sys.argv[3]

	outfile = open( output_file_name, "w" )

	infile2 = open( input_file_2, "r" )
	inlines2 = infile2.readlines()
	infile2.close()
	inlines2.sort()

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
	inlines1.sort()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			file2_seq = return_seq_for_id( this_id )

			if (file2_seq != this_seq):
				output_line = this_id + "\r\n"
				outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

