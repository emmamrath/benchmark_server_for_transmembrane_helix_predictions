#!/usr/bin/env python

# python list_ids_in_both_files.py opm_feb2012.all_opm.opm_TMsubunit opm_feb2012.alpha_polytopic.opm_seq opm_feb2012.alpha_polytopic.opm_seq.and.opm_TMsubunit

# python list_ids_in_file1_not_in_file2.py [input-file-1] [input-file-2] [output-file]

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

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			this_id = get_id( inline )
			found_id = look_for_id( this_id )
			if (found_id == 1):
				output_line = this_id + "\r\n"
				outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

