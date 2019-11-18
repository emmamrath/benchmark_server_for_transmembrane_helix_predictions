#!/usr/bin/env python

# python ../remove_HHHHHH.py opm_feb2012.list2c_from_pdb_seq.pdb_seq opm_feb2012.list2d_from_pdb_seq.pdb_seq

# python rremove_HHHHHH.py [input-file [output-file]

import sys
import re


# global variables
save_pdbid = []
save_id = []
save_id_in_not_found_array = []


######################################################
def remove_HHHHHH( this_seq ):

	new_seq = this_seq
	len_this_seq = len(this_seq)
	if ( len_this_seq >= 6 ):
		last_6_residues = new_seq[ (len_this_seq-6) : len_this_seq ]
		if ( last_6_residues == 'HHHHHH' ):
			new_seq = ''
			in_HHHHHH = 1
			for i in range( (len_this_seq-1), -1, -1 ):
				this_aa = this_seq[i:(i+1)]
				if (in_HHHHHH == 1):
					if (this_aa == 'H'):
						do_nothing = 1
					else:
						in_HHHHHH = 0
				if (in_HHHHHH == 0):
					new_seq = this_aa + new_seq
	return new_seq

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
			new_seq = remove_HHHHHH( this_seq )
			output_line = this_id + ':::::' + new_seq + ':::::' + "\r\n"
			outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

