#!/usr/bin/env python

# python convert_1line_field3_to_fasta.py opm_feb2012_beta_barrel.pdb_seq.missing_blanks opm_feb2012_beta_barrel.pdb_seq.missing_blanks.fasta
# python convert_1line_field3_to_fasta.py [input-file] [output-file]

import sys
import re


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

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_fasta_id = bits[1]
			this_seq = bits[2]
			output_line = this_fasta_id + "\r\n"
			outfile.write( output_line )
			output_line = this_seq + "\r\n"
			outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

