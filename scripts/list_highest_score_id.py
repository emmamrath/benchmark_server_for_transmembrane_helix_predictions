#!/usr/bin/env python

# python ../list_highest_score_id.py list11_not_in_mpstruc_pass1_pass2_pass3.needle.235_mpstruc_ids.similarity_scores list11_not_in_mpstruc_pass1_pass2_pass3.needle.235_mpstruc_ids.similarity_scores.highest_score

# python list_highest_score_id.py [input-file] [output-file]

# Example of input file:
# - : 1a11_A,1a91_A,1afo_A
# 1h68_A : 28.5,53.3,47.4,
# 2ei4_A : 45.7,43.2,71.1,
# 1py6_A : 34.4,47.6,54.4,

import sys
import re

# global variables
save_id = []


######################################################
def main():

	global outfile
	global output_file_name
	input_file = sys.argv[1]
	output_file_name = sys.argv[2]

	outfile = open( output_file_name, "w" )

	infile = open( input_file, "r" )
	inlines = infile.readlines()
	infile.close()

	inline_header = inlines[0]
	bits = inline_header.rsplit( ' : ' )
	inline_header = bits[1]
	header_ids = inline_header.rsplit( ',' )
	for this_id in header_ids:
		this_id = this_id.strip()
		save_id.append( this_id )

	j = 0
	for inline in inlines:
		inline = inline.strip()
		if ((inline != '') and (j != 0)):
			bits = inline.rsplit( ' : ' )
			this_id = bits[0]
			this_id = this_id.strip()
			list_of_scores = bits[1]
			scores = list_of_scores.rsplit( ',' )
			i = 0
			highest_score = 0
			highest_id = ''
			for this_score in scores:
				this_score = this_score.strip()
				if (this_score > highest_score):
					highest_score = this_score
					highest_id = save_id[i]
				i = i + 1
			if (highest_id == ''):
				highest_id = '-'
			output_line = this_id + ':::::' + highest_id + ':::::' + highest_score + ':::::' + "\r\n"
			outfile.write( output_line )
		j = j + 1
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

