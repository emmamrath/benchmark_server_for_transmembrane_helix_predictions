#!/usr/bin/env python

# INCOMPLETE PROGRAM

# python compare_seqs_for_position_equality.py 3bz1.opm_seq 3bz1.pdb_seq

# python compare_seqs_for_position_equality.py [input-file-1] [input-file-2] [output-file]

import sys
import re


# global variables
save_id = []
save_seq = []
found_seq = ''


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
def get_seq( inline ):

	bits = inline.rsplit( ':::::' )
	this_seq = bits[1]
	return this_seq

######################################################
def look_for_id( this_id ):

	global found_seq
	found_id = 0
	i = 0
	for this_save_id in save_id:
		if (this_save_id == this_id):
			found_id = 1
			found_seq = save_seq[i]
		i = i + 1
	return found_id

######################################################
def compare_seqs( seq1, seq2 ):

	if (len(seq1) < len(seq2)):
		for i in range( len(seq1), len(seq2) ):
			seq1 = seq1 + '_'
	else:
		if (len(seq2) < len(seq1)):
			for i in range( len(seq2), len(seq1) ):
				seq2 = seq2 + '_'
	seqs_have_same_positions = 1
	for i in range( 0, len(seq1) ):
		j = i + 1
		if ((seq1[i:j] == '_') or (seq2[i:j] == '_')):
			dont_compare_this_char = 1
		else:
			if (seq1[i:j] != seq2[i:j]):
				seqs_have_same_positions = 0
	return seqs_have_same_positions

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

	for inline in inlines2:
		inline = inline.strip()
		if (inline != ''):
			this_id = get_id( inline )
			this_seq = get_seq( inline )
			if (this_seq != '-'):
				save_id.append( this_id )
				save_seq.append( this_seq )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			this_id = get_id( inline )
			found_id = look_for_id( this_id )
			if (found_id == 1):
				this_seq = get_seq( inline )
				seqs_have_same_positions = compare_seqs( this_seq, found_seq )
				if (seqs_have_same_positions == 0):
					output_line = 'SEQS NOT THE SAME FOR ' + this_id + " : \r\n"
					outfile.write( output_line )
					output_line = '     file 1 : ' + this_seq + "\r\n"
					outfile.write( output_line )
					output_line = '     file 2 : ' + found_seq + "\r\n"
					outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

