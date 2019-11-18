#!/usr/bin/env python

# python sort_file_contents_by_file2.py alpha_polytopic_bitopic.pdb_seq alpha_polytopic_bitopic.earliest_release_date.sort1 alpha_polytopic_bitopic.pdb_seq.sort1

# python sort_file_contents_by_file2.py [input-file-1] [input-file-2] [output-file]

import sys
import re


# global variables
save_id = []
save_record = []


######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	output_file_name = sys.argv[3]

	outfile = open( output_file_name, "w" )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			this_id = inline[0:6]
			save_id.append( this_id )
			save_record.append( inline )

	infile2 = open( input_file_2, "r" )
	inlines2 = infile2.readlines()
	infile2.close()

	for inline in inlines2:
		inline = inline.strip()
		if (inline != ''):
			this_id = inline[0:6]
			this_save_record = ''
			i = 0
			for this_save_id in save_id:
				if (this_save_id == this_id):
					this_save_record = save_record[i]
				i = i + 1
			if (this_save_record == ''):
				print 'Did not find record for', this_id
			else:
				output_line = this_save_record + "\r\n"
				outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

