#!/usr/bin/env python

# python ../list_file1_field2_ids_having_file2_pdbid.py in_list11_not_in_mpstruc.same_seq_ids mpstruc_category.pdbid in_list11_not_in_mpstruc.mpstruc_pdbid

# python list_file1_field2_ids_having_file2_pdbid.py [input-file-1] [input-file-2] [output-file]

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
			bits = inline.rsplit( ':::::' )
			this_pdbid = bits[0]
			save_pdbid.append( this_pdbid )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_list_of_ids = bits[1]
			list_of_ids = this_list_of_ids.rsplit( ',' )
			found_pdbid = 0
			for this_id_in_list in list_of_ids:
				this_pdbid_in_list = this_id_in_list[0:4]
				for this_save_pdbid in save_pdbid:
					if (this_pdbid_in_list == this_save_pdbid):
						if (found_pdbid == 0):
							found_pdbid = 1
							output_line = this_id + ':::::' + this_save_pdbid + ':::::' + "\r\n"
							outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

