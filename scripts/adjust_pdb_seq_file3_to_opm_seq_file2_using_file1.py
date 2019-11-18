#!/usr/bin/env python

# python ../adjust_pdb_seq_file3_to_opm_seq_file2_using_file1.py alpha_polytopic_bitopic_peptide.opm_seq_not_equal_pdb_seq.unequal_start_posn alpha_polytopic_bitopic_peptide.opm_seq alpha_polytopic_bitopic_peptide.pdb_seq alpha_polytopic_bitopic_peptide.opm_seq_not_equal_pdb_seq.renum_pdb_seq

# python adjust_pdb_seq_file3_to_opm_seq_file2_using_file1.py [input-file-1-list-renumbering-info] [input-file-2-opm-seq] [input-file-3-pdb-seq] [output-file-pdb-seq-renumbered-according-to-file1]

import sys
import re


# global variables
save_opm_id = []
save_opm_seq = []
save_pdb_id = []
save_pdb_seq = []


######################################################
def return_opm_seq_for_id( this_id ):

	found_seq = ''
	i = 0
	for this_save_id in save_opm_id:
		if (this_save_id == this_id):
			found_seq = save_opm_seq[i]
		i = i + 1
	return found_seq

######################################################
def return_pdb_seq_for_id( this_id ):

	found_seq = ''
	i = 0
	for this_save_id in save_pdb_id:
		if (this_save_id == this_id):
			found_seq = save_pdb_seq[i]
		i = i + 1
	return found_seq

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	input_file_3 = sys.argv[3]
	output_file_name = sys.argv[4]

	outfile = open( output_file_name, "w" )

	infile2 = open( input_file_2, "r" )
	inlines2 = infile2.readlines()
	infile2.close()
	for inline in inlines2:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			save_opm_id.append( this_id )
			save_opm_seq.append( this_seq )

	infile3 = open( input_file_3, "r" )
	inlines3 = infile3.readlines()
	infile3.close()
	for inline in inlines3:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			save_pdb_id.append( this_id )
			save_pdb_seq.append( this_seq )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_pdb_start = bits[1]
			this_opm_start = bits[2]
			this_pdb_seq = return_pdb_seq_for_id( this_id )
			this_opm_seq = return_opm_seq_for_id( this_id )

			new_pdb_seq = this_pdb_seq

			if (this_pdb_start > this_opm_start):
				i = int(this_pdb_start) - 1
				new_pdb_seq = this_pdb_seq[i:]

			if (this_pdb_start < this_opm_start):
				new_pdb_seq = ''
				for i in range( 1, int(this_opm_start) ):
					new_pdb_seq = '_' + new_pdb_seq
				new_pdb_seq = new_pdb_seq + this_pdb_seq

			output_line = this_id + ':::::' + new_pdb_seq + ':::::' + "\r\n"
			outfile.write( output_line )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

