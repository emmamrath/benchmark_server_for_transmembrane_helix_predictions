#!/usr/bin/env python

# python ../output_file1_id_with_field2_match_file2_pdbid.py in_list11_not_in_mpstruc.found_same_seq_mpstruc.mpstruc_pdbid mpstruc_feb2012.pdbid.mpstruc_category list11_in_mpstruc_pass1.mpstruc_category

# python output_file1_id_with_field2_match_file2_pdbid.py [input-file-1] [input-file-2] [output-file]

import sys
import re


# global variables
save_pdbid = []
save_field2 = []


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
			bits = inline.rsplit( ':::::' )
			this_pdbid = bits[0]
			this_field2 = bits[1]
			save_pdbid.append( this_pdbid )
			save_field2.append( this_field2 )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_pdbid = this_id[0:4]
			have_done_output = 0
			i = 0
			for this_save_pdbid in save_pdbid:
				if (have_done_output == 0):
					if (this_save_pdbid == this_pdbid):
						this_save_field2 = save_field2[i]
						output_line = this_id + ':::::' + this_save_field2 + ':::::' + "\r\n"
						outfile.write( output_line )
						have_done_output = 1
				i = i + 1
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

