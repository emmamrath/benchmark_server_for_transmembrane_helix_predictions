#!/usr/bin/env python

# python expand_pdbid_file_contents_id_list.py opm_feb2012_beta_barrel.pdbid.pdb_method_resolution opm_feb2012_beta_barrel.id opm_feb2012_beta_barrel.pdb_method_resolution.expanded

# python expand_pdbid_file_contents_id_list.py [input-file-1] [input-file-2] [output-file]

import sys
import re


# global variables
save_id = []


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
			this_id = inline[0:6]
			save_id.append( this_id )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			this_pdbid = inline[0:4]
			for this_save_id in save_id:
				this_save_pdbid = this_save_id[0:4]
				if (this_save_pdbid == this_pdbid):
					bits = inline.split( ':::::' )
					this_record = bits[1]
					output_line = this_save_id + ':::::' + this_record + ':::::' + "\r\n"
					outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

