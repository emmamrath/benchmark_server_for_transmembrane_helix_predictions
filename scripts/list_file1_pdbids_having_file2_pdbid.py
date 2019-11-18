#!/usr/bin/env python

# python list_file1_pdbids_having_file2_pdbid.py opm_feb2012.all_opm_having_pdb_struc2D.pdb_seq possible_missing_chains_in_pdb_struc2D

# python list_file1_pdbids_having_file2_pdbid.py [input-file-1] [input-file-2] [output-file]

import sys
import re


# global variables
save_pdbid = []
save_id = []
output_id = []
output_pdbid = []


######################################################
def get_pdbid( inline ):

	this_pdbid = inline[0:4]
	return this_pdbid

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
def get_chain( inline ):

	this_chain = inline[5:6]
	return this_chain

######################################################
def is_it_lowercase( in_str ):

	it_is_lower = 0
	if (in_str == in_str.lower()):
		it_is_lower = 1
	return it_is_lower

######################################################
def look_for_pdbid( this_pdbid ):

	found_pdbid = 0
	for this_save_pdbid in save_pdbid:
		if (this_save_pdbid == this_pdbid):
			found_pdbid = 1
	return found_pdbid

######################################################
def look_for_id( this_id ):

	found_id = 0
	for this_save_id in save_id:
		if (this_save_id == this_id):
			found_id = 1
	return found_id

######################################################
def look_for_output_id( this_id ):

	found_id = 0
	for this_output_id in output_id:
		if (this_output_id == this_id):
			found_id = 1
	return found_id

######################################################
def look_for_output_pdbid( this_pdbid ):

	found_pdbid = 0
	for this_output_pdbid in output_pdbid:
		if (this_output_pdbid == this_pdbid):
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
			this_id = get_id( inline )
			save_id.append( this_id )
			this_pdbid = get_pdbid( inline )
			save_pdbid.append( this_pdbid )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			this_pdbid = get_pdbid( inline )
			i = 0
			for this_save_pdbid in save_pdbid:
				if (this_save_pdbid == this_pdbid):
					found_it = look_for_output_pdbid( this_save_pdbid )
					if (found_it == 0):
						output_line = this_save_pdbid + "\r\n"
						outfile.write( output_line )
						output_pdbid.append( this_save_pdbid )
				i = i + 1
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

