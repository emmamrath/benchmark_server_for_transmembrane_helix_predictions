#!/usr/bin/env python

# python align_file1_seq_to_file2.py 3bz1.pdb_seq 3bz1.opm_seq 3bz1.pdb_seq_aligned

# python align_file1_seq_to_file2.py [input-file-1] [input-file-2] [output-file]

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
def compare_seqs_wrong( seq1, seq2 ):

	#if (len(seq1) < len(seq2)):
	#	for i in range( len(seq1), len(seq2) ):
	#		seq1 = seq1 + '_'
	#else:
	#	if (len(seq2) < len(seq1)):
	#		for i in range( len(seq2), len(seq1) ):
	#			seq2 = seq2 + '_'

	compare_seq2 = ''
	in_seq2_blanks = 1
	for i in range( 0, len(seq2) ):
		if ((seq2[i] == '_') and (in_seq2_blanks == 1)):
			do_nothing = 1
		else:
			in_seq2_blanks = 0
			compare_seq2 = compare_seq2 + seq2[i]

	seqs_have_same_positions = 1
	for i in range( 0, len(compare_seq2) ):
		j = i + 1
		if ((seq1[i:j] == '_') or (compare_seq2[i:j] == '_')):
			dont_compare_this_char = 1
		else:
			if (seq1[i:j] != compare_seq2[i:j]):
				seqs_have_same_positions = 0
	return seqs_have_same_positions

######################################################
def compare_seqs( seq1, seq2 ):

	if (len(seq2) > len(seq1)):
		seqs_have_same_positions = 0
	else:
		seqs_have_same_positions = 1
		for i in range( 0, len(seq2) ):
			j = i + 1
			if ((seq1[i:j] == '_') or (seq2[i:j] == '_')):
				dont_compare_this_char = 1
			else:
				if (seq1[i:j] != seq2[i:j]):
					seqs_have_same_positions = 0
	return seqs_have_same_positions

######################################################
def align_seq1_to_seq2( seq1, seq2, this_id ):

	new_seq1 = seq1

	found_start_posn_seq2 = 0
	reached_end_of_seq2 = 0
	seq2_posn = 0
	start_posn_seq2 = 0
	while ((found_start_posn_seq2 == 0) and (reached_end_of_seq2 == 0)):
		if ( seq2[seq2_posn:(seq2_posn+1)] != '_' ):
			found_start_posn_seq2 = 1
			start_posn_seq2 = seq2_posn
		else:
			seq2_posn = seq2_posn + 1
			if (seq2_posn >= len(seq2)):
				reached_end_of_seq2 = 1

	if (found_start_posn_seq2 == 1):
		found_aligning_posn = 0
		reached_end_of_seq1 = 0
		try_aligning_seq1_posn = 0
		seq1_align_posn = 0
		seq2_align_posn = 0
		while ((found_aligning_posn == 0) and (reached_end_of_seq1 == 0)):
			this_posn_is_aligning_nicely = 1
			reached_end_of_either_seq = 0
			seq1_posn = try_aligning_seq1_posn
			seq2_posn = start_posn_seq2
			while ((this_posn_is_aligning_nicely == 1) and (reached_end_of_either_seq == 0)):
				if ( seq1[seq1_posn:(seq1_posn+1)] == seq2[seq2_posn:(seq2_posn+1)] ):
					seq1_posn = seq1_posn + 1
					seq2_posn = seq2_posn + 1
					if ((seq1_posn >= len(seq1)) or (seq2_posn >= len(seq2))):
						reached_end_of_either_seq = 1
				else:
					this_posn_is_aligning_nicely = 0
			if (this_posn_is_aligning_nicely == 1):
				found_aligning_posn = 1
				seq1_align_posn = try_aligning_seq1_posn
				seq2_align_posn = start_posn_seq2
				new_seq1 = ''
				for i in range( 0, (seq2_align_posn - seq1_align_posn)):
					new_seq1 = new_seq1 + '_'
				new_seq1 = new_seq1 + seq1
				# Seq2 came from reading residue numbers in a PDB file.
				# If the PDB file does not start at residue number 1,
				# then seq2 contains '_' from residue number 1.
				# Seq1 came from PDB file's FASTA sequence.
				# Sometimes the PDB file's FASTA sequence starts after residue number 1 too.
				# In that case, the output sequence will have '_' 
				# because we don't have the beginning residues in either the FASTA sequence or the PDB file residues.
				# However, there are some opposite cases where the PDB file starts at residue 1,
				# and yet the PDB FASTA sequence has residues before that.
				# This is probably an error in the PDB files.
				# I need to flag them, because the OPM TMH segments will be numbered according to the PDB file.
				# However, I would like to submit the entire FASTA sequence to TMH prediction programs, especially HMMs.
				# So I will flag these cases, and renumber residues in the OPM TMH segments the same as the FASTA sequence.
				#
				# if (seq2_align_posn < seq1_align_posn):
				#	print 'NUMBERING PROBLEM :', this_id
				#
				# I won't flag them after all. There are too many of them to fix.
				# Instead, truncate the output sequence to match the PDB file and OPM TMH segments numbering.
				if (seq2_align_posn < seq1_align_posn):
					new_seq1 = new_seq1[ (seq1_align_posn - seq2_align_posn): ]
			else:
				try_aligning_seq1_posn = try_aligning_seq1_posn + 1
				if (try_aligning_seq1_posn >= len(seq1)):
					reached_end_of_seq1 = 1

	return new_seq1

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
			save_id.append( this_id )
			save_seq.append( this_seq )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			print 'processing', this_id
			this_id = get_id( inline )
			found_id = look_for_id( this_id )
			if (found_id == 1):
				this_seq = get_seq( inline )
				seqs_have_same_positions = compare_seqs( this_seq, found_seq )
				if (seqs_have_same_positions == 1):
					output_line = inline + "\r\n"
					outfile.write( output_line )
				else:
					aligned_seq = align_seq1_to_seq2( this_seq, found_seq, this_id )
					output_line = this_id + ':::::' + aligned_seq + ":::::\r\n"
					outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

