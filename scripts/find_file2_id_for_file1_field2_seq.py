#!/usr/bin/env python

# python ../find_file2_id_for_file1_field2_seq.py in_list10_not_in_pdbtm.pdb_seq pdbtm_pdbids_from_pdbtm_having_pdb_seq.pdb_seq seq_in_list10_seq_in_pdbtm

# python find_file2_id_for_file1_field2_seq.py [input-file-1] [input-file-2] [output-file]

import sys
import re


# global variables
save_id = []
save_seq = []


######################################################
def look_for_seq( this_seq ):

	found_id = ''
	i = 0
	for this_save_seq in save_seq:
		if (this_save_seq == this_seq):
			found_id = save_id[i]
		i = i + 1
	return found_id

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

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			found_id = look_for_seq( this_seq )
			if (found_id == ''):
				output_line = this_id + ':::::' + 'AN ID FOR SEQ NOT FOUND' + ':::::' + "\r\n"
			else:
				output_line = this_id + ':::::' + found_id + ':::::' + "\r\n"
			outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

