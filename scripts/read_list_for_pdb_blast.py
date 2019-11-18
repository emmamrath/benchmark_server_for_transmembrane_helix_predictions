#!/usr/bin/env python

# python read_list_for_pdb_blast.py alpha_polytopic_bitopic.pdb_seq alpha_polytopic_bitopic.earliest_homologue_release_date

# python read_list_for_pdb_blast.py [input-list-of-ids-and-sequences] [output-file-for-dates]

# This program reads a list of PDBID+chain and calls PDB BLAST.

import sys
import os
import re
import SOAPpy


######################################################
# global variables


######################################################
def main():

	global outfile
	input_list_file = sys.argv[1]
	output_file_name = sys.argv[2]

	outfile = open( output_file_name, "w" )
	outfile.close()

	infile = open( input_list_file, "r" )
	inlines = infile.readlines()
	infile.close()

	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			this_seq = this_seq.replace( '_', '' )
			print 'processing', this_id
			this_command = 'python do_one_pdb_blast.py ' + this_id + ' ' + this_seq + ' ' + output_file_name
			os.system( this_command )


if __name__ == "__main__":
	# Someone is launching this directly
	main()

