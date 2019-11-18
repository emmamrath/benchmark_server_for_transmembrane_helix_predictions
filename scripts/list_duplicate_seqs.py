#!/usr/bin/env python

# python list_duplicate_seqs.py opm_feb2012.opm_category_list8_alpha_polytopic_bitopic_peptide.pdb_seq opm_feb2012.opm_category_list8_alpha_polytopic_bitopic_peptide.pdb_seq.dupl_seqs

# python list_duplicate_seqs.py [input-file] [output-file]

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
def look_for_seq( this_seq ):

	found_it = 0
	for this_save_seq in save_seq:
		if (this_save_seq == this_seq):
			found_it = 1
	return found_it

######################################################
def look_for_dupl_seq( this_seq ):

	found_it = 0
	for this_save_seq in save_dupl_seq:
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
				save_seq.append( seq )
			else:
				found_it = look_for_dupl_seq( seq )
				if (found_it == 0):
					save_dupl_seq.append( seq )
					output_line =  seq + "\r\n"
					outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

