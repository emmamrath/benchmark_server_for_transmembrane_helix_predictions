#!/usr/bin/env python

# python ../read_pdb_seq_for_list_pdbid.py one_mpstruc_category.pdbid one_mpstruc_category.pdb_seq

# python read_pdb_seq_for_list_pdbid.py [input-file] [output-file]

import sys
import re
import SOAPpy


# global variables
save_seq = []


######################################################
def get_pdb_seqs( this_pdbid, server ):

	global outfile
	global output_file_name

	formatted_line = ''

	print 'processing', this_pdbid

	soappy_chains_output = server.getChains ( this_pdbid );

	for this_chain in soappy_chains_output:

		this_pdbid_and_chain = this_pdbid + '_' + this_chain

		soappy_output = server.getSequenceForStructureAndChain ( this_pdbid, this_chain );

		this_seq = soappy_output

		if (this_seq == 'N/A'):
			print 'PROBLEM :', this_pdbid, this_chain, 'returned', this_seq
			this_seq = '-'
		else:

			found_seq = 0
			for this_save_seq in save_seq:
				if (this_save_seq == this_seq):
					found_seq = 1

			if (found_seq == 0):
				save_seq.append( this_seq )
				output_line = this_pdbid_and_chain + ':::::' + this_seq + ':::::' + "\r\n"
				outfile.write( output_line )

	outfile.close()
	outfile = open( output_file_name, "a" )

	return

######################################################
def main():

	global outfile
	global output_file_name
	input_file = sys.argv[1]
	output_file_name = sys.argv[2]

	outfile = open( output_file_name, "w" )

	server = SOAPpy.SOAPProxy("http://www.pdb.org/pdb/services/pdbws")

	infile = open( input_file, "r" )
	inlines = infile.readlines()
	infile.close()

	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			this_pdbid = inline[0:4]
			get_pdb_seqs( this_pdbid, server )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

