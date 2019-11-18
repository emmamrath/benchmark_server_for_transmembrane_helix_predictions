#!/usr/bin/env python

# python read_one_pdb_struc2D.py 1bl8_A 1bl8_A.pdb_struc2D

# python read_one_pdb_struc2D.py [input-pdbid-and-chain] [output-append-file-for-struc2D]

# This program makes a call to PDB to get the 2D structure for a PDBID+chain.

import sys
import os
import re
import SOAPpy


######################################################
# global variables

######################################################
def get_pdb_struc2D( input_id, server ):

	formatted_line = ''

	this_pdbid = input_id[0:4]
	this_pdbid_uc = this_pdbid.upper()
	this_chain = input_id[5:6]
	this_pdbid_and_chain = this_pdbid + '_' + this_chain

	soappy_output = server.getKabschSander ( this_pdbid_uc, this_chain );

	struc2D = soappy_output

	struc2D = struc2D.replace( ' ', '_' )

	if (struc2D == 'N/A'):
		print 'PROBLEM :', this_pdbid_uc, this_chain, 'returned', struc2D
		struc2D = '-'

	formatted_line = this_pdbid_and_chain + ':::::' + struc2D + ':::::'

	return formatted_line

######################################################
def write_output_struc2D( output_line, output_file_name ):

	outfile = open( output_file_name, "a" ) # open output file for writing in append mode

	output_line = output_line + "\r\n"
	outfile.write( output_line )

	outfile.close()

	return


######################################################
def main():

	global outfile
	input_id = sys.argv[1]
	output_file_name = sys.argv[2]

	server = SOAPpy.SOAPProxy("http://www.pdb.org/pdb/services/pdbws")

	output_line = get_pdb_struc2D( input_id, server )

	write_output_struc2D( output_line, output_file_name )


if __name__ == "__main__":
	# Someone is launching this directly
	main()

