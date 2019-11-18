#!/usr/bin/env python

# python sort_file.py opm_feb2012.list1_from_pdb_seq.pdb_seq opm_feb2012.sorted_list1_from_pdb_seq.pdb_seq

# python sort_file.py [input-file] [output-file]

import sys
import re


# global variables


######################################################
def get_id( inline ):

	formatted_line = ''

	#this_pdbid = inline[0:4]
	#this_chain = inline[5:6]
	#this_pdbid_and_chain = this_pdbid + '_' + this_chain
	this_pdbid_and_chain = inline[0:6]

	formatted_line = this_pdbid_and_chain

	return formatted_line

######################################################
def main():

	global outfile
	input_file = sys.argv[1]
	output_file_name = sys.argv[2]

	outfile = open( output_file_name, "w" )

	infile = open( input_file, "r" )
	inlines = infile.readlines()
	infile.close()

	inlines.sort()

	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			output_line = inline + "\r\n"
			outfile.write( output_line )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

