#!/usr/bin/env python

# python ../../../read_list_opm_pdb_for_stride.py list_opm_alpha_polytopic_files.txt opm_feb2012.alpha_polytopic.stride_seq opm_feb2012.stride_struc2D

# python read_list_opm_pdb_for_stride.py [input-list-of-pdb-files] [output-file-for-segments] [output-file-for-sequences] [output-file-for-membrane-near]

# This program reads a list of PDB files from OPM where the membrane position has been defined by DUM records.
# It calls STRIDE for each file to get secondary structure, and outputs the sequence and secondary structure to one output file each.

import sys
import os


######################################################
# global variables


######################################################
def main():

	global outfile
	input_list_file = sys.argv[1]
	output_file_name_seq = sys.argv[2]
	output_file_name_struc2D = sys.argv[3]

	outfile_seq = open( output_file_name_seq, "w" )
	outfile_struc2D = open( output_file_name_struc2D, "w" )
	outfile_seq.close()
	outfile_struc2D.close()

	infile = open( input_list_file, "r" )
	inlines = infile.readlines()
	infile.close()

	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			this_pdb_file = inline
			this_command = 'python ../call_stride_for_one_opm_pdb.py ' + this_pdb_file + ' ' + output_file_name_seq + ' ' + output_file_name_struc2D
			os.system( this_command )


if __name__ == "__main__":
	# Someone is launching this directly
	main()

