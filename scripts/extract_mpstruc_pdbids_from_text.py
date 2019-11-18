#!/usr/bin/env python

# python extract_mpstruc_pdbids_from_text.py mpstruc_expanded_list.txt mpstruc_pdbid_list.txt

# python extract_mpstruc_pdbids_from_text.py [input-file] [output-file]

import sys
import re


# global variables

######################################################
def get_id( inline ):

	this_id = inline[0:4]
	return this_id

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
			if (len(inline) == 4):
				output_line = inline + "\r\n"
				outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

