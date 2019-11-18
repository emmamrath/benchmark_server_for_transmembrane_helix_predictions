#!/usr/bin/env python

# python assign_category.py opm_feb2012.alpha_bitopic.opm_seq 'OPM-TRANSMEMBRANE-HELIXL;;;;;OPM-BITOPIC;;;;;' opm_feb2012.alpha_bitopic.opm_category

# python list_ids_in_file1_not_in_file2.py [input-file] [input-category] [output-file]

import sys
import re


# global variables
save_id = []


######################################################
def get_id( inline ):

	this_id = inline[0:6]
	#this_pdbid = inline[0:4]
	#this_pdbid = this_pdbid.lower()
	#this_chain = inline[5:6]
	#this_chain = this_chain.upper()
	#this_id = this_pdbid + '_' + this_chain
	return this_id

######################################################
def main():

	global outfile
	input_file = sys.argv[1]
	input_category = sys.argv[2]
	output_file_name = sys.argv[3]

	outfile = open( output_file_name, "w" )

	infile = open( input_file, "r" )
	inlines = infile.readlines()
	infile.close()

	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			this_id = get_id( inline )
			output_line = this_id + ':::::' + input_category + ':::::' + "\r\n"
			outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

