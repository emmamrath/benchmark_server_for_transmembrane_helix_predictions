#!/usr/bin/env python

# python parse_scampi_single_topologies.py alpha_polytopic_bitopic_noblanks__scampi.cbr.su.se__topologies.txt alpha_polytopic_bitopic.scampi

# python parse_scampi_single_topologies.py [input-file] [output-file]

import sys
import re


# input file :
# >2b6o_A
# iiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMoMMMMMMMMMMMMMMMMMMMMMi
# iiiiiiiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooooooooooo
# oooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiMMMMMMMMMMMMMMMMMMMMMooo
# oooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiii
# iiiiiiiiiiiiiiiiiiiiiii
# 
# >2zt9_G
# ooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiii


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

	in_a_seq = 0
	curr_id = ''
	curr_seq = ''
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			if (inline[0:1] == '>'):

				if (in_a_seq == 1):
					if (curr_seq == ''):
						curr_seq = '-'
					output_line = curr_id + ':::::' + curr_seq + ':::::' + "\r\n"
					outfile.write( output_line )

				curr_id = inline[1:7]
				curr_seq = ''
				in_a_seq = 1
			else:
				curr_seq = curr_seq + inline

	if (in_a_seq == 1):
		if (curr_seq == ''):
			curr_seq = '-'
		output_line = curr_id + ':::::' + curr_seq + ':::::' + "\r\n"
		outfile.write( output_line )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

