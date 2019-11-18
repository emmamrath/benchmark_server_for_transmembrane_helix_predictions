#!/usr/bin/env python

# python read_pdb_enzymeclassification_for_list.py opm_feb2012.all_opm.opm_seq opm_feb2012.all_opm.pdb_ecnum

# python read_pdb_enzymeclassification_for_list.py [input-file] [output-file]

import sys
import re
import SOAPpy


# global variables


######################################################
def get_pdb_info( inline, server ):

	formatted_line = ''

	this_pdbid = inline[0:4]
	this_pdbid_uc = this_pdbid.upper()
	this_chain = inline[5:6]
	this_pdbid_and_chain = this_pdbid + '_' + this_chain

	print 'processing ', this_pdbid_and_chain

	soappy_output = server.getEcNumsForStructures ( this_pdbid_uc );

	pdbinfo = soappy_output
	print soappy_output
	# <<class 'SOAPpy.Types.typedArrayType'> getEcNumsForStructuresReturn at 151907660>: []

	if (pdbinfo == 'N/A'):
		print 'PROBLEM :', inline, 'returned', pdbinfo
		pdbinfo = '-'

	formatted_line = this_pdbid_and_chain + ':::::' + pdbinfo + ':::::'

	return formatted_line

######################################################
def main():

	global outfile
	input_file = sys.argv[1]
	output_file_name = sys.argv[2]

	server = SOAPpy.SOAPProxy("http://www.pdb.org/pdb/services/pdbws")

	outfile = open( output_file_name, "w" )

	infile = open( input_file, "r" )
	inlines = infile.readlines()
	infile.close()

	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			output_line = get_pdb_info( inline, server )
			output_line = output_line + "\r\n"
			outfile.write( output_line )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

