#!/usr/bin/env python

# python remove_duplicate_ids.py opm_feb2012_in_pdb.pdb_struc2D.got_duplicates opm_feb2012_in_pdb.pdb_struc2D

# python remove_duplicate_ids.py [input-file] [output-file]

import sys
import re


# global variables
save_id = []


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
def look_for_id( this_id ):

	found_id = 0
	for this_save_id in save_id:
		if (this_save_id == this_id):
			found_id = 1
	return found_id

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
			this_id = get_id( inline )
			found_id = look_for_id( this_id )
			if (found_id == 0):
				output_line = inline + "\r\n"
				outfile.write( output_line )
				save_id.append( this_id )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

