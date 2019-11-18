#!/usr/bin/env python

# python list_ids_from_pdbselect.py recent.pdb_select25.txt recent.pdb_select25.id

# python list_ids_from_pdbselect.py [input-file] [output-file]

import sys
import re

#    25  1KDXB    28 -1.00  0.00   NMR    28    28      1     21      0 ' creb'
#    25  1JXCA    68 -1.00  0.00   NMR    68    68      0      9      9 ' putative trypsin inhibitor atti-2'
#    25  1K4UP    32 -1.00  0.00   NMR    32    32      0     12      0 ' phagocyte nadph oxidase subunit p47phox'



# global variables


######################################################
def main():

	global outfile
	input_file = sys.argv[1]
	output_file_name = sys.argv[2]

	outfile = open( output_file_name, "w" )

	infile = open( input_file, "r" )
	inlines = infile.readlines()
	infile.close()

	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			inline = inline.replace( '  ', ' ' )
			inline = inline.replace( '  ', ' ' )
			inline = inline.replace( '  ', ' ' )
			inline = inline.replace( '  ', ' ' )
			bits = inline.rsplit( ' ' )
			this_pdbid_chain = bits[1]
			this_pdbid = this_pdbid_chain[0:4]
			this_pdbid = this_pdbid.lower()
			this_chain = this_pdbid_chain[4:5]
			this_id = this_pdbid + '_' + this_chain
			output_line =  this_id + "\r\n"
			outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

