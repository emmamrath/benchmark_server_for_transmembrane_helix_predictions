#!/usr/bin/env python

# python list_unique_pdbids.py opm_feb2012.unique_real_seq_opm_in_pdb.id opm_feb2012.unique_real_seq_opm_in_pdb.pdbid

# python list_unique_pdbid.py [input-file] [output-file]

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
save_pdbid = []

######################################################
def look_for_pdbid( this_pdbid ):

	found_it = 0
	for this_save_pdbid in save_pdbid:
		if (this_save_pdbid == this_pdbid):
			found_it = 1
	return found_it

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
			pdbid = inline[0:4]
			found_it = look_for_pdbid( pdbid )
			if (found_it == 0):
				save_pdbid.append( pdbid )
				output_line = pdbid + "\r\n"
				outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

