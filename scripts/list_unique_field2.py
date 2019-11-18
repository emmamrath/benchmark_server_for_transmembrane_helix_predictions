#!/usr/bin/env python

# python list_unique_field2.py opm_feb2012.unique_real_seq_opm_in_pdb.id opm_feb2012.unique_real_seq_opm_in_pdb.pdbid

# python list_unique_field2.py [input-file]

import sys
import re

# input file :
# 1gzm_A
# 1h2s_A
# 1h2s_B
# 1h68_A

# output file :
# 1gzm
# 1h2s
# 1h68

# global variables
save_field = []

######################################################
def process_field( this_field ):

	if ((this_field != '') and (this_field != '-')):
		found_it = 0
		for this_save_field in save_field:
			if (this_save_field == this_field):
				found_it = 1
		if (found_it == 0):
			save_field.append( this_field )
	return

######################################################
def main():

	global outfile
	input_file = sys.argv[1]
	#output_file_name = sys.argv[2]

	infile = open( input_file, "r" )
	inlines = infile.readlines()
	infile.close()

	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			field2 = bits[1]
			bits2 = field2.rsplit( ';;;;;' )
			for field2_subsection in bits2:
				bits3 = field2_subsection.rsplit( ',,,,,' )
				field2_name = bits3[0]
				process_field( field2_name )

	save_field.sort()

	#outfile = open( output_file_name, "w" )
	for this_field in save_field:
		output_line = this_field + "\r\n"
		print this_field
		#outfile.write( output_line )
	#outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

