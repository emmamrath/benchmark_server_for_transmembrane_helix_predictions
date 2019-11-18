#!/usr/bin/env python

# python ../list_ids.py list11_alpha_polytopic_bitopic_peptide.pdb_seq list11_alpha_polytopic_bitopic_peptide.id

# python list_ids.py [input-file] [output-file]

import sys
import re

# 1fjk_A:::::MDKVQYLTRSAIRRASTIEMPQQARQNLQNLFINFCLILIFLLLICIIVMLL:::::
# 1grm_A:::::VGALAVVVWLWLWLWX:::::
# 1grm_B:::::VGALAVVVWLWLWLWX:::::
# 1grm_A:::::VGALAVVVWLWLWLWX:::::


# global variables
save_seq = []
save_dupl_seq = []

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
			this_id = bits[0]
			output_line =  this_id + "\r\n"
			outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

