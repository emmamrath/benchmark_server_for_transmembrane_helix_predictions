#!/usr/bin/env python

# python parse_topcons_single_topologies.py alpha_polytopic_bitopic_noblanks__single.topcons.net__topologies.txt alpha_polytopic_bitopic.topcons

# python parse_topcons_single_topologies.py [input-file] [output-file]

import sys
import re


# input file :
# >3arc_K
# iiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooo
# 
# >2e74_B
# iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMoooooooooooooooooooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooooooooooo
# 
# >1qle_D
# iiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooo


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

	curr_id = ''
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			if (inline[0:1] == '>'):
				curr_id = inline[1:7]
			else:
				output_line = curr_id + ':::::' + inline + ':::::' + "\r\n"
				outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

