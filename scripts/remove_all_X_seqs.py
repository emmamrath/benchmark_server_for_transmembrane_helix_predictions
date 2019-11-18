#!/usr/bin/env python

# python ../remove_all_X_seqs.py one_mpstruc_category.id.pdb_seq one2_mpstruc_category.id.pdb_seq

# python remove_all_X_seqs.py [input-file] [output-file]

# Remove X from file 1.
# Remove X at the end of a sequence.
# Turn X into _ in the middle of a sequence.
# Remove X at the beginning of a sequence. Compare it to file2. If X is the first position, then replace it with _.
#	If the next residue is really the first residue, then just remove the X.

import sys
import re


# global variables
save_id = []
save_seq = []
found_seq = ''


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
			this_seq = bits[1]
			found_non_X = 0
			for i in range( 0, len(this_seq) ):
				j = i + 1
				this_char = this_seq[i:j]
				if (this_char != 'X'):
					found_non_X = 1
			if (found_non_X == 1):
				output_line = inline + "\r\n"
				outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

