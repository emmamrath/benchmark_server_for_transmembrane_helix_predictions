#!/usr/bin/env python

# python ../list_ids_in_file1_not_in_pdbid_file2.py list11_alpha_polytopic_bitopic_peptide.pdb_seq mpstruc_pdbid_list.txt in_list11_not_in_mpstruc

# python list_ids_in_file1_not_in_pdbid_file2.py.py [input-file-1] [input-file-2] [output-file]

import sys
import re


# global variables
save_pdbid = []
save_pdbid_in_not_found_array = []


######################################################
def get_id( inline ):

	this_id = inline[0:6]
	return this_id

######################################################
def get_pdbid( inline ):

	this_pdbid = inline[0:4]
	return this_pdbid

######################################################
def look_for_pdbid( this_pdbid ):

	found_pdbid = 0
	for this_save_pdbid in save_pdbid:
		if (this_save_pdbid == this_pdbid):
			found_pdbid = 1
	return found_pdbid

######################################################
def look_for_pdbid_in_not_found_array( this_pdbid ):

	found_pdbid = 0
	for this_save_pdbid in save_pdbid_in_not_found_array:
		if (this_save_pdbid == this_pdbid):
			found_pdbid = 1
	return found_pdbid

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	output_file_name = sys.argv[3]

	outfile = open( output_file_name, "w" )

	infile2 = open( input_file_2, "r" )
	inlines2 = infile2.readlines()
	infile2.close()

	for inline in inlines2:
		inline = inline.strip()
		if (inline != ''):
			this_pdbid = get_pdbid( inline )
			found_pdbid = look_for_pdbid( this_pdbid )
			if (found_pdbid == 0):
				save_pdbid.append( this_pdbid )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			this_id = get_id( inline )
			this_pdbid = get_pdbid( inline )
			found_pdbid = look_for_pdbid( this_pdbid )
			if (found_pdbid == 0):
				output_line = this_id + "\r\n"
				outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

