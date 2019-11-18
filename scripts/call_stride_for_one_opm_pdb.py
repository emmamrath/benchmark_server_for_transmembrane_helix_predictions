#!/usr/bin/env python

# python call_stride_for_one_opm_pdb.py 1r3j.pdb 1r3j.stride_seq 1r3j.stride_struc2D

# python call_stride_for_one_opm_pdb.py [input-pdb-file] [output-append-file-for-stride-seq] [output-append-file-for-stride-struc2D]

# This program reads a PDB file from OPM where the membrane position has been defined by DUM records.
# It calls STRIDE to get the secondary structure

import sys
import os


######################################################
# global variables
resn_list       = ['ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','GLY','PRO','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']
resn_list_1char = ['R',  'H',  'K',  'D',  'E',  'S',  'T',  'N',  'Q',  'C',  'G',  'P',  'A',  'V',  'I',  'L',  'M',  'F',  'Y',  'W']
chain_array = []
resi_array = []
resn_array = []
struc2D_array = []


######################################################
def is_residue( resn ):

	this_is_residue = 0
	for resn_in_array in resn_list:
		if (resn_in_array == resn):
			this_is_residue = 1
	return this_is_residue


######################################################
def resn_3char_to_1char( resn_3char ):

	resn_1char = 'X'
	for i in range( 0, len(resn_list) ):
		resn_in_array = resn_list[i]
		if (resn_in_array == resn_3char):
			resn_1char = resn_list_1char[i]
	return resn_1char


######################################################
def read_stride_file( input_file ):

	# *)  IMPORTANT NOTE: if the protein chain	identifier is '	' (space), it
	#     will	be substituted by '-' (dash) everywhere	in the STRIDE output.
	#     The same is true  for  command  line	 parameters  involving	chain
	#     identifiers where you have to specify '-' instead of	' '.

	# **) One-letter secondary	structure code is nearly the same as used  in
	#     DSSP	[2] (see Frishman and Argos [1]	for details):
	# 
	#    H	    Alpha helix
	#    G	    3-10 helix
	#    I	    PI-helix
	#    E	    Extended conformation
	#    B or	b   Isolated bridge
	#    T	    Turn
	#    C	    Coil (none of the above)

	# SEQ  1    WIQITLQP                                              8          1R9U
	# STR       TTTTT                                                            1R9U
	# REM                                                                        1R9U
	# REM                                                                        1R9U
	# REM                                                                        1R9U
	# LOC  TurnIV       TRP     2 A      ILE      6 A                            1R9U
	# LOC  TurnIV       ILE     3 A      THR      7 A                            1R9U
	# REM                                                                        1R9U
	# REM  --------------- Detailed secondary structure assignment-------------  1R9U
	# REM                                                                        1R9U
	# REM  |---Residue---|    |--Structure--|   |-Phi-|   |-Psi-|  |-Area-|      1R9U
	# ASG  TRP A    2    1    T          Turn    360.00    -34.89     232.7      1R9U
	# ASG  ILE A    3    2    T          Turn    -70.95    -25.12      95.0      1R9U
	# ASG  GLN A    4    3    T          Turn    -77.25    360.00     191.3      1R9U
	# ASG  ILE A    6    4    T          Turn    360.00    -33.84     110.2      1R9U
	# ASG  THR A    7    5    T          Turn    -68.30    360.00     128.5      1R9U
	# ASG  LEU A    9    6    C          Coil    360.00    360.00     160.0      1R9U
	# ASG  GLN A   12    7    C          Coil    360.00    360.00     219.2      1R9U
	# ASG  PRO A   16    8    C          Coil    360.00    360.00     220.4      1R9U

	# ASG  ASN C   50   50    E        Strand   -155.32    173.37       1.4      3VGA
	# ASG  ILE C   51   51    E        Strand   -145.50    138.17       0.8      3VGA
	# ASG  ASN C   52   52    T          Turn    -87.60    114.76      40.2      3VGA
	# ASG  PRO C  53A   53    T          Turn    -72.07    -16.38       0.2      3VGA

	infile = open( input_file, "r" )
	inlines = infile.readlines()
	infile.close()

	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			this_record_type = inline[0:3]
			if (this_record_type == 'ASG'):
				this_chain = inline[9:10]
				this_resn_3char = inline[5:8]
				this_resn_1char = resn_3char_to_1char(this_resn_3char)
				this_resi = inline[10:15]
				this_resi_5 = inline[14:15]
				is_char = 0
				try:
					x = int(this_resi_5)
					is_char = 0
				except ValueError:
					is_char = 1
				if (is_char == 1):
					this_resi = inline[10:14]
				this_resi = this_resi.strip()
				this_resi = int(this_resi)
				this_struc2d = inline[24:25]

				chain_array.append( this_chain )
				resi_array.append( this_resi )
				resn_array.append( this_resn_1char )
				struc2D_array.append( this_struc2d )
	return


######################################################
def write_output_seq_and_struc2D( this_pdbid, output_file_name_seq, output_file_name_struc2D ):

	outfile_seq = open( output_file_name_seq, "a" ) # open output file for writing in append mode
	outfile_struc2D = open( output_file_name_struc2D, "a" ) # open output file for writing in append mode

	prev_chain = ''
	prev_resi = 0
	seq_string = ''
	struc2D_string = ''

	seen_chain_array = []

	for i in range( 0, len(chain_array) ):

		this_chain = chain_array[i]
		this_resi = resi_array[i]
		this_resn = resn_array[i]
		this_struc2D = struc2D_array[i]

		if (this_chain != prev_chain):

			# at the end of reading a chain, save the entire chain sequence
			seen_this_chain_before = 0
			if (prev_chain != ''):

				# output the previous chain
				this_id = this_pdbid + '_' + prev_chain
				output_line_seq = this_id + ':::::' + seq_string + ':::::' + "\r\n"
				output_line_struc2D = this_id + ':::::' + struc2D_string + ':::::' + "\r\n"
				outfile_seq.write( output_line_seq )
				outfile_struc2D.write( output_line_struc2D )

				seen_this_chain_before = 0
				for this_seen_chain in seen_chain_array:
					if (this_seen_chain == this_chain):
						seen_this_chain_before = 1

				seq_string = ''
				struc2D_string = ''
				prev_resi = 0

			# if we've already seen this chain before, then this second appearance is probably an amino acid ligand entered as 'ATOM  '
			if (seen_this_chain_before == 0):

				seen_chain_array.append( this_chain )
				prev_resi_plus_1 = prev_resi + 1
				while (prev_resi_plus_1 < this_resi):
					seq_string = seq_string + '_'
					struc2D_string = struc2D_string + '_'
					prev_resi_plus_1 = prev_resi_plus_1 + 1

				seq_string = seq_string + this_resn
				struc2D_string = struc2D_string + this_struc2D

		else:

			prev_resi_plus_1 = prev_resi + 1
			while (prev_resi_plus_1 < this_resi):
				seq_string = seq_string + '_'
				struc2D_string = struc2D_string + '_'
				prev_resi_plus_1 = prev_resi_plus_1 + 1

			seq_string = seq_string + this_resn
			struc2D_string = struc2D_string + this_struc2D

		prev_resi = this_resi
		prev_chain = this_chain

	if (seq_string != ''):

		# output the last chain
		this_id = this_pdbid + '_' + prev_chain
		output_line_seq = this_id + ':::::' + seq_string + ':::::' + "\r\n"
		output_line_struc2D = this_id + ':::::' + struc2D_string + ':::::' + "\r\n"
		outfile_seq.write( output_line_seq )
		outfile_struc2D.write( output_line_struc2D )

	outfile_seq.close()
	outfile_struc2D.close()

	return


######################################################
def main():

	global outfile
	input_pdb_file = sys.argv[1]
	output_file_name_seq = sys.argv[2]
	output_file_name_struc2D = sys.argv[3]

	this_pdbid = input_pdb_file[0:4]
	temp_stride_output_file = 'stride.out'

	this_command = '/home/emma/Emmas_files_on_linuxbox/software_to_install/stride/stride/stride ' + input_pdb_file + ' -f' + temp_stride_output_file
	os.system( this_command )

	read_stride_file( temp_stride_output_file )

	#print 'chain_array =', chain_array
	#print 'resi_array =', resi_array
	#print 'resn_array =', resn_array
	#print 'struc2D_array =', struc2D_array

	write_output_seq_and_struc2D( this_pdbid, output_file_name_seq, output_file_name_struc2D )


if __name__ == "__main__":
	# Someone is launching this directly
	main()

