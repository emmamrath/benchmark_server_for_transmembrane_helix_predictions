#!/usr/bin/env python

# python compare_similarity_of_lists_of_2_lists.py 'alpha_polytopic_bitopic.pdb_seq.sort1.water.similarity_scores.*.binary_matrix.nonhomologous_list' 'sort2/alpha_polytopic_bitopic.pdb_seq.sort2.water.similarity_scores.*.binary_matrix.nonhomologous_list'

# python compare_similarity_of_lists_of_2_lists.py [input-file-1] [input-file-2]

import sys
import re
import os

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
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]

	i = 20
	while (i <= 100):
		file1 = input_file_1.replace( '*', str(i) )
		file2 = input_file_2.replace( '*', str(i) )
		
		this_command = 'python compare_similarity_of_2_lists.py ' + str(i) + ' ' + file1 + ' ' + file2
		os.system( this_command )

		i = i + 5



if __name__ == "__main__":
	# Someone is launching this directly
	main()

