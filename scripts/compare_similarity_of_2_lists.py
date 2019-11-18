#!/usr/bin/env python

# python compare_similarity_of_2_lists.py 20 alpha_polytopic_bitopic.pdb_seq.sort1.water.similarity_scores.20.binary_matrix.nonhomologous_list sort2/alpha_polytopic_bitopic.pdb_seq.sort2.water.similarity_scores.20.binary_matrix.nonhomologous_list

# compare_similarity_of_2_lists.py [input-file-1] [input-file-2]

import sys
import re

# 1fjk_A:::::MDKVQYLTRSAIRRASTIEMPQQARQNLQNLFINFCLILIFLLLICIIVMLL:::::
# 1grm_A:::::VGALAVVVWLWLWLWX:::::
# 1grm_B:::::VGALAVVVWLWLWLWX:::::
# 1grm_A:::::VGALAVVVWLWLWLWX:::::


# global variables
save_file1_id = []
save_file2_id = []


######################################################
def main():

	global outfile
	input_id = sys.argv[1]
	input_file_1 = sys.argv[2]
	input_file_2 = sys.argv[3]

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	num_ids = 0
	num_same_ids = 0

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			this_id = inline[0:6]
			save_file1_id.append( this_id )
			num_ids = num_ids + 1

	infile2 = open( input_file_2, "r" )
	inlines2 = infile2.readlines()
	infile2.close()

	for inline in inlines2:
		inline = inline.strip()
		if (inline != ''):
			this_id = inline[0:6]
			save_file2_id.append( this_id )

	for this_save_file1_id in save_file1_id:
		for this_save_file2_id in save_file2_id:
			if (this_save_file1_id == this_save_file2_id):
				num_same_ids = num_same_ids + 1

	#print num_same_ids, num_ids
	percent_same = int((float(num_same_ids) / num_ids) * 100)
	print 'For', input_id, '%', 'list, percent same is', percent_same, '%'

	return

if __name__ == "__main__":
	# Someone is launching this directly
	main()

