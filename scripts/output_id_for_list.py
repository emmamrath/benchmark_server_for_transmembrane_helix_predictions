#!/usr/bin/env python

# python output_id_for_list.py opm_feb2012.pdb_seq opm_feb2012.id

# python output_id_for_list.py [input-file] [output-file]

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

	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			output_line = get_id( inline )
			output_line = output_line + "\r\n"
			outfile.write( output_line )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

