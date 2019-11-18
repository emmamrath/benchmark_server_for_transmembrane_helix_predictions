#!/usr/bin/env python

# python list_unique_seqs.py mpstruc_expanded_list.txt mpstruc_pdbid_list.txt

# python list_unique_seqs.py [input-file] [output-file]

import sys
import re

# 1fjk_A:::::MDKVQYLTRSAIRRASTIEMPQQARQNLQNLFINFCLILIFLLLICIIVMLL:::::
# 1grm_A:::::VGALAVVVWLWLWLWX:::::
# 1grm_B:::::VGALAVVVWLWLWLWX:::::
# 1grm_A:::::VGALAVVVWLWLWLWX:::::


# global variables
save_inline = []
save_seq = []

######################################################
def look_for_seq( this_seq ):

	found_it = 0
	for this_save_seq in save_seq:
		if (this_save_seq == this_seq):
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
			bits = inline.rsplit( ':::::' )
			seq = bits[1]
			found_it = look_for_seq( seq )
			if (found_it == 0):
				save_inline.append( inline )
				save_seq.append( seq )
				output_line = inline + "\r\n"
				outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

