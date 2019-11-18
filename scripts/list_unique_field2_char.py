#!/usr/bin/env python

# python ../list_unique_field2_char.py alpha_polytopic_bitopic_peptide.pdb_struc2D

# python list_unique_field2_char.py [input-file]

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
save_char = []

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
			for i in range( 0, len(field2) ):
				j = i + 1
				this_char = field2[i:j]
				is_saved_already = 0
				for this_save_char in save_char:
					if (this_save_char == this_char):
						is_saved_already = 1
				if (is_saved_already == 0):
					save_char.append( this_char )

	save_char.sort()

	for this_save_char in save_char:
		print this_save_char


if __name__ == "__main__":
	# Someone is launching this directly
	main()

