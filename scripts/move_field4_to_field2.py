#!/usr/bin/env python

# python move_field4_to_field2.py ../pdbtm_feb2012.all_pdbtm.TMH in_pdbtm_not_in_pdb_opm.txt

# python move_field4_to_field2.py [input-file-1] [output-file]

import sys
import re


# global variables


######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	output_file_name = sys.argv[2]

	outfile = open( output_file_name, "w" )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_field4 = bits[3]
			output_line = this_id + ':::::' + this_field4 + ':::::' + "\r\n"
			outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

