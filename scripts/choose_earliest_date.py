#!/usr/bin/env python

# python choose_earliest_date.py alpha_polytopic_bitopic.earliest_homologue_release_date alpha_polytopic_bitopic.release_date alpha_polytopic_bitopic.earliest_release_date

# python choose_earliest_date.py [input-file-1] [input-file-2] [output-file]

import sys
import re


# global variables
save_id = []
save_date = []


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
			this_id = bits[0]
			this_date = bits[1]
			save_id.append( this_id )
			save_date.append( this_date )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_date = bits[1]
			i = 0
			found_it = 0
			this_save_date = '9999-99-99'
			for this_save_id in save_id:
				if (this_save_id == this_id):
					found_it = 1
					this_save_date = save_date[i]
				i = i + 1

			output_line = this_id + ':::::' + this_date + ':::::' + "\r\n"
			if (this_save_date < this_date):
				output_line = this_id + ':::::' + this_save_date + ':::::' + "\r\n"				
			outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

