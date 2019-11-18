#!/usr/bin/env python

# python list_file1_pdbids_not_in_file2_id.py OPM_beta_barrel_proteins.pdbid opm_feb2012_from_zip.beta_barrel.id OPM_beta_barrel_not_got_file

# python list_file1_pdbids_not_in_file2_id.py [input-file-1] [input-file-2] [output-file]

import sys
import re


# global variables
save_pdbid = []


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
			this_pdbid = inline[0:4]
			save_pdbid.append( this_pdbid )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			this_id = inline[0:6]
			this_pdbid = inline[0:4]
			found_it = 0
			for this_save_pdbid in save_pdbid:
				if (this_save_pdbid == this_pdbid):
					found_it = 1
			if (found_it == 0):
				output_line = this_pdbid + "\r\n"
				outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

