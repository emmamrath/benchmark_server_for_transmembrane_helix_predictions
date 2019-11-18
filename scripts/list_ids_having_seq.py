#!/usr/bin/env python

# python list_ids_having_seq.py opm_feb2012.opm_category_list8_alpha_polytopic_bitopic_peptide.pdb_seq.dupl_seqs opm_feb2012.opm_category_list8_alpha_polytopic_bitopic_peptide.pdb_seq opm_feb2012.opm_category_list8_alpha_polytopic_bitopic_peptide.ids_for_dupl_seqs

# python list_ids_having_seq.py [input-file-1] [input-file-2] [output-file]

import sys
import re

# 1fjk_A:::::MDKVQYLTRSAIRRASTIEMPQQARQNLQNLFINFCLILIFLLLICIIVMLL:::::
# 1grm_A:::::VGALAVVVWLWLWLWX:::::
# 1grm_B:::::VGALAVVVWLWLWLWX:::::
# 1grm_A:::::VGALAVVVWLWLWLWX:::::


# global variables
save_id = []
save_seq = []
output_seq = []
output_ids = []

######################################################
def add_id_to_seq( this_seq, this_id ):

	i = 0
	found_it = 0
	for this_output_seq in output_seq:
		if (this_output_seq == this_seq):
			this_output_ids = output_ids[i]
			this_output_ids = this_output_ids + this_id + ','
			output_ids[i] = this_output_ids
			found_it = 1
		i = i + 1
	if (found_it == 0):
		output_seq.append( this_seq )
		output_ids.append( this_id + ',' )
	return

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	output_file_name = sys.argv[3]

	infile2 = open( input_file_2, "r" )
	inlines2 = infile2.readlines()
	infile2.close()
	for inline in inlines2:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			save_id.append( this_id )
			save_seq.append( this_seq )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_seq = bits[0]
			i = 0
			for this_save_seq in save_seq:
				if (this_save_seq == this_seq):
					this_save_id = save_id[i]
					add_id_to_seq( this_seq, this_save_id )
				i = i + 1

	outfile = open( output_file_name, "w" )
	i = 0
	for this_seq in output_seq:
		this_ids = output_ids[i]
		output_line =  this_seq + ':::::' + this_ids + ':::::' + "\r\n"
		outfile.write( output_line )
		i = i + 1
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

