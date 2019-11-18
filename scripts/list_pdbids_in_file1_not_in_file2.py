#!/usr/bin/env python

# python ../list_pdbids_in_file1_not_in_file2.py opm_feb2012.all_opm.pdb_seq ../PDBTM/pdbtm_feb2012.all_TM_prot.TMH pdbid_in_pdb_opm_not_in_pdbtm.txt
# python ../list_pdbids_in_file1_not_in_file2.py ../PDBTM/pdbtm_feb2012.all_TM_prot.TMH opm_feb2012.all_opm.pdb_seq pdbid_in_pdbtm_not_in_pdb_opm.txt

# python list_pdbids_in_file1_not_in_file2.py [input-file-1] [input-file-2] [output-file]

import sys
import re


# global variables
save_pdbid = []
save_pdbid_in_not_found_array = []


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
			this_pdbid = get_pdbid( inline )
			found_pdbid = look_for_pdbid( this_pdbid )
			if (found_pdbid == 0):
				found_pdbid_in_not_found_array = look_for_pdbid_in_not_found_array( this_pdbid )
				if (found_pdbid_in_not_found_array == 0):
					# output_line = 'PDBID NOT FOUND : ' + this_pdbid + "\r\n"
					output_line = this_pdbid + "\r\n"
					outfile.write( output_line )
					save_pdbid_in_not_found_array.append( this_pdbid )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

