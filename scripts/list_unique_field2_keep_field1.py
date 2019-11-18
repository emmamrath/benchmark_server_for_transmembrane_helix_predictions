#!/usr/bin/env python

# python ../list_unique_field2_keep_field1.py mpstruc_category.pdbid one_mpstruc_category.pdbid

# python list_unique_field2_keep_field1.py [input-file-1] [output-file]

import sys
import re


# global variables
save_pdbid = []


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

	prev_mpstruc_category = ''
	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_mpstruc_category = bits[1]
			if (prev_mpstruc_category != this_mpstruc_category):
				output_line = this_id + ':::::' + this_mpstruc_category + ':::::' + "\r\n"
				outfile.write( output_line )
				prev_mpstruc_category = this_mpstruc_category
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

