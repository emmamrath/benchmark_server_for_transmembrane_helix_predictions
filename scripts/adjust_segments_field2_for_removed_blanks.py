#!/usr/bin/env python

# python adjust_segments_field2_for_removed_blanks.py test1.fasta.split4 test1.leading_blanks_removed test1.fasta.split4.adjusted_for_removed_blanks

# python adjust_segments_field2_for_removed_blanks.py [input-file] [input-file-2] [output-file]

import sys
import re


# global variables
save_id = []
save_num_removed = []


######################################################
def find_num_removed( this_id ):
	this_num_removed = -1
	i = 0
	for this_save_id in save_id:
		if (this_save_id == this_id):
			this_num_removed = save_num_removed[i]
		i = i + 1
	return this_num_removed

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	output_file_name = sys.argv[3]

	infile2 = open( input_file_2, "r" )
	inlines2 = infile2.readlines()
	infile2.close()

	for inline in inlines2:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_num_removed = bits[1]
			save_id.append( this_id )
			save_num_removed.append( this_num_removed )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	outfile = open( output_file_name, "w" )

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_segments = bits[1]
			this_num_removed = find_num_removed( this_id )
			this_num_removed = int(this_num_removed)
			if (this_num_removed == -1):
				print this_id, ': could not find a record for how many blanks were removed.'
			else:
				output_line = this_id + ':::::' + this_segments + ':::::' + "\r\n"
				if (this_num_removed > 0):
					new_segments = ''
					if ((this_segments != '-') and (this_segments != '-,')):
						bits2 = this_segments.rsplit( ',' )
						for this_segment in bits2:
							if (this_segment != ''):
								bits3 = this_segment.rsplit( '-' )
								this_from = bits3[0]
								this_to = bits3[1]
								this_from = int(this_from) + this_num_removed
								this_to = int(this_to) + this_num_removed
								new_segment = str(this_from) + '-' + str(this_to) + ','
								new_segments = new_segments + new_segment
					if (new_segments == ''):
						new_segments = '-'
					output_line = this_id + ':::::' + new_segments + ':::::' + "\r\n"
				outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

