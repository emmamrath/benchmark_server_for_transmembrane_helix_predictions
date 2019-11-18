#!/usr/bin/env python

# python parse_tmhmm2_long_web_page.py soluble.tmhmm2_results.txt soluble.tmhmm2

# python parse_tmhmm2_long_web_page.py [input-file] [output-file]

import sys
import re


# input file :
# # 2r4r_A Length: 365
# # 2r4r_A Number of predicted TMHs:  7
# # 2r4r_A Exp number of AAs in TMHs: 150.1817
# # 2r4r_A Exp number, first 60 AAs:  22.97602
# # 2r4r_A Total prob of N-in:        0.00399
# # 2r4r_A POSSIBLE N-term signal sequence
# 2r4r_A	TMHMM2.0	outside	     1    35
# 2r4r_A	TMHMM2.0	TMhelix	    36    58
# 2r4r_A	TMHMM2.0	inside	    59    69
# 2r4r_A	TMHMM2.0	TMhelix	    70    92
# 2r4r_A	TMHMM2.0	outside	    93   106
# 2r4r_A	TMHMM2.0	TMhelix	   107   129
# 2r4r_A	TMHMM2.0	inside	   130   149
# 2r4r_A	TMHMM2.0	TMhelix	   150   169
# 2r4r_A	TMHMM2.0	outside	   170   200
# 2r4r_A	TMHMM2.0	TMhelix	   201   223
# 2r4r_A	TMHMM2.0	inside	   224   274
# 2r4r_A	TMHMM2.0	TMhelix	   275   297
# 2r4r_A	TMHMM2.0	outside	   298   306
# 2r4r_A	TMHMM2.0	TMhelix	   307   326
# 2r4r_A	TMHMM2.0	inside	   327   365

# # 1vkk_A Length: 154
# # 1vkk_A Number of predicted TMHs:  0
# # 1vkk_A Exp number of AAs in TMHs: 0.00666
# # 1vkk_A Exp number, first 60 AAs:  0
# # 1vkk_A Total prob of N-in:        0.22994
# 1vkk_A	TMHMM2.0	outside	     1   154


# global variables

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	output_file_name = sys.argv[2]

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	outfile = open( output_file_name, "w" )

	prev_id = ''
	this_id = ''
	this_seq = ''
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):

			inline_is_data_line = 0
			if ((inline[0:1] != '_') and (inline[4:5] == '_')):
				inline_is_data_line = 1

			if (inline_is_data_line == 1):

				this_id = inline[0:6]

				if ((this_id != prev_id) and (prev_id != '')):
					output_line = prev_id + ':::::' + this_seq + ':::::' + "\r\n"
					outfile.write( output_line )

				if (this_id != prev_id):
					this_seq = ''	
					this_char = '_'

				if ( inline.find('TMhelix') != -1):
					this_char = 'M'
				if ( inline.find('inside') != -1):
					this_char = 'i'
				if ( inline.find('outside') != -1):
					this_char = 'o'

				inline = inline.replace( "\t", ' ' )
				inline = inline.replace( '  ', ' ' )
				inline = inline.replace( '  ', ' ' )
				inline = inline.replace( '  ', ' ' )
				inline = inline.replace( '  ', ' ' )
				inline = inline.replace( '  ', ' ' )
				inline = inline.replace( '  ', ' ' )
				bits = inline.rsplit( ' ' )

				this_from = int(bits[3])
				this_to = int(bits[4])
				this_to_plus_1 = this_to + 1

				for j in range( this_from, this_to_plus_1 ):
					this_seq = this_seq + this_char

				prev_id = this_id

	if (this_id != ''):
		output_line = this_id + ':::::' + this_seq + ':::::' + "\r\n"
		outfile.write( output_line )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

