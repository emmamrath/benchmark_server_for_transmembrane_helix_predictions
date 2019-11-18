#!/usr/bin/env python

# python list_file1_ids_having_file2_id.py opm_feb2012.all_opm.opm_seq possible_wrong_upper_chains_in_pdb_struc2D need_to_redo_pdb_struc2D_for_lowercase

# python list_file1_ids_having_file2_id.py [input-file-1] [input-file-2] [output-file]

import sys
import re


# global variables
save_pdbid = []
save_id = []
output_id = []
output_pdbid = []


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
			this_id = inline[0:6]
			save_id.append( this_id )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			this_id = inline[0:6]
			i = 0
			for this_save_id in save_id:
				if (this_save_id == this_id):
					output_line = this_id + "\r\n"
					outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

