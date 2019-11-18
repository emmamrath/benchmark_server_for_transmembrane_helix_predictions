#!/usr/bin/env python

# python ../list_not_equal_seqs_for_non_missing_residues.py alpha_polytopic_bitopic_peptide.opm_seq.sorted alpha_polytopic_bitopic_peptide.renum_pdb_seq alpha_polytopic_bitopic_peptide.opm_seq_not_equal_renum_pdb_seq_for_non_missing_residues.id

# python list_not_equal_seqs_for_non_missing_residues.py [input-file-1] [input-file-2] [output-file]

import sys
import re


# global variables
save_id = []
save_seq = []


######################################################
def return_seq_for_id( this_id ):

	found_seq = ''
	i = 0
	for this_save_id in save_id:
		if (this_save_id == this_id):
			found_seq = save_seq[i]
		i = i + 1
	return found_seq

######################################################
def judge_if_similar( seq1, seq2 ):

	seqs_are_similar = 1

	len_seq1 = len(seq1)
	len_seq2 = len(seq2)
	smallest_len = len_seq1
	if (len_seq2 < len_seq1):
		smallest_len = len_seq2

	for i in range( 0, smallest_len ):
		j = i + 1
		char_seq1 = seq1[i:j]
		char_seq2 = seq2[i:j]
		if ((char_seq1 != '_') and (char_seq2 != '_') and (char_seq1 != 'X') and (char_seq2 != 'X')):
			if (char_seq1 != char_seq2):
				seqs_are_similar = 0

	return seqs_are_similar

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	output_file_name = sys.argv[3]

	outfile = open( output_file_name, "w" )

	infile2 = open( input_file_2, "r" )
	inlines2 = infile2.readlines()
	infile2.close()
	inlines2.sort()

	for inline in inlines2:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			save_id.append( this_id )
			save_seq.append( this_seq )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()
	inlines1.sort()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			file2_seq = return_seq_for_id( this_id )

			seqs_are_similar = judge_if_similar( this_seq, file2_seq )

			if (seqs_are_similar == 0):
				output_line = this_id + "\r\n"
				outfile.write( output_line )
				#print this_id
				#print this_seq
				#print file2_seq
				#print
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

