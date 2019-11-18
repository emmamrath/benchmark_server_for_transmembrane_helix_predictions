#!/usr/bin/env python

# python ../list_same_seq_ids_for_list.py in_list11_not_in_mpstruc.pdb_seq opm_feb2012.list4b.pdb_seq in_list11_not_in_mpstruc.same_seq_ids

# python list_same_seq_ids_for_list.py [input-file-1] [input-file-2] [output-file]

import sys
import re


# global variables
save_id = []
save_seq = []


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
			list_of_ids = ''
			i = 0
			for this_save_seq in save_seq:
				if (this_save_seq == this_seq):
					this_save_id = save_id[i]
					list_of_ids = list_of_ids + this_save_id + ','
				i = i + 1
			output_line = this_id + ':::::' + list_of_ids + ':::::' + "\r\n"
			outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

