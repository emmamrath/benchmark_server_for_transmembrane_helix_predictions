#!/usr/bin/env python

# python ../compare_opm_pdbtm_topology.py alpha_polytopic_bitopic.opm_adj_tm_segments_topology alpha_polytopic_bitopic.pdbtm_segments_topology

# python compare_opm_pdbtm_topology.py [input-file]

import sys
import re


# global variables
save_id = []
save_file2 = []


######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]

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
			save_file2.append( this_seq )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			found_it = 0
			i = 0
			#this_save_seq = ''
			for this_save_id in save_id:
				if (this_save_id == this_id):
					found_it = 1
					this_save_seq = save_file2[i]
				i = i + 1
			if (found_it == 0):
				print 'Did not find a match in file2 for file1', this_id
			else:
				mismatch = 0
				for i in range( 0, len(this_seq) ):
					j = i + 1
					if ((this_seq[i:j] == 'i') and (this_save_seq[i:j] == 'o')):
						mismatch = 1
					if ((this_seq[i:j] == 'o') and (this_save_seq[i:j] == 'i')):
						mismatch = 1
				if (mismatch == 1):
					print 'There is a topology mismatch for', this_id


if __name__ == "__main__":
	# Someone is launching this directly
	main()

