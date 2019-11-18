#!/usr/bin/env python

# python read_list_pdb_struc2D.py test1.txt test1.pdb_struc2D

# python read_list_pdb_struc2D.py [input-list-of-ids] [output-file-for-struc2D]

# This program reads a list of PDBID+chain and calls PDB to get the struc2D.

import sys
import os
import re
import SOAPpy


######################################################
# global variables


######################################################
def get_id( inline ):

	this_id = inline[0:6]
	#this_pdbid = inline[0:4]
	#this_pdbid = this_pdbid.lower()
	#this_chain = inline[5:6]
	#this_chain = this_chain.upper()
	#this_id = this_pdbid + '_' + this_chain
	return this_id


######################################################
def main():

	global outfile
	input_list_file = sys.argv[1]
	output_file_name = sys.argv[2]

	outfile = open( output_file_name, "w" )

	infile = open( input_list_file, "r" )
	inlines = infile.readlines()
	infile.close()

	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			this_id = get_id( inline )
			print 'processing', this_id
			this_command = 'python ../read_one_pdb_struc2D.py ' + this_id + ' ' + output_file_name
			os.system( this_command )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

