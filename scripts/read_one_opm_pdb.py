#!/usr/bin/env python

# python read_one_opm_pdb.py 1r3j.pdb alpha_polytopic.noseq.opm_feb2012 opm_feb2012.alpha_polytopic.membrane_posn opm_feb2012.alpha_polytopic.seq opm_feb2012.alpha_polytopic.opm_membrane_near

# python read_one_opm_pdb.py [input-pdb-file] [output-append-file-for-segments] [output-append-file-for-sequences] [output-append-file-for-membrane-near]

# This program reads a PDB file from OPM where the membrane position has been defined by DUM records.
# It outputs the segments as being M (membrane), I (inside the cell or on the inside side of the membrane)
# or O (outside the cell or on the outside side of the membrane).
# This program does not judge the secondary structure of the segments,
# which could be alpha helices, beta sheets or coils.

import sys


######################################################
# global variables
outfile = ''
resn_list       = ['ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','GLY','PRO','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']
resn_list_1char = ['R',  'H',  'K',  'D',  'E',  'S',  'T',  'N',  'Q',  'C',  'G',  'P',  'A',  'V',  'I',  'L',  'M',  'F',  'Y',  'W']
chain_array = []
resi_array = []
resn_array = []
z_array = []
posn_array = []
membrane_near_array = []
DUM_N_z = 0
DUM_O_z = 0
got_DUM_N = 0
got_DUM_O = 0
membrane_segments_segments = []
membrane_segments_chain = []
inside_segments_segments = []
inside_segments_chain = []
outside_segments_segments = []
outside_segments_chain = []
all_segments_segments = []
all_segments_chain = []
all_segments_membrane_near = []


# ATOM      1  N   ALA A  23     -16.378   6.310 -19.289  1.00181.62           N  
# ATOM      2  CA  ALA A  23     -15.241   5.381 -19.555  1.00181.62           C  
# ATOM      3  C   ALA A  23     -15.333   4.149 -18.663  1.00181.62           C  
# ATOM      4  O   ALA A  23     -16.317   3.409 -18.705  1.00181.62           O  
# ATOM      5  CB  ALA A  23     -13.914   6.100 -19.317  1.00 74.09           C  
# ATOM      6  N   LEU A  24     -14.296   3.936 -17.858  1.00163.39           N  
# ATOM      7  CA  LEU A  24     -14.240   2.806 -16.940  1.00163.39           C  
# ATOM      8  C   LEU A  24     -13.111   3.009 -15.938  1.00163.39           C  
# ATOM      9  O   LEU A  24     -13.348   3.278 -14.765  1.00163.39           O  
# ATOM     10  CB  LEU A  24     -14.010   1.503 -17.709  1.00126.45           C  
# ATOM     11  CG  LEU A  24     -13.979   0.232 -16.856  1.00126.45           C  
# ATOM     12  CD1 LEU A  24     -15.243  -0.576 -17.100  1.00126.45           C  
# ATOM     13  CD2 LEU A  24     -12.744  -0.586 -17.192  1.00126.45           C  
# HETATM 2146  N   DUM  2146     -26.000   0.000 -16.200
# HETATM 2146  O   DUM  2146     -26.000   0.000  16.200
# HETATM 2222  N   DUM  2222     -24.000 -10.000 -16.200
# HETATM 2222  O   DUM  2222     -24.000 -10.000  16.200
# HETATM 2223  N   DUM  2223     -24.000  -8.000 -16.200
# HETATM 2223  O   DUM  2223     -24.000  -8.000  16.200
# HETATM 2224  N   DUM  2224     -24.000  -6.000 -16.200
# HETATM 2224  O   DUM  2224     -24.000  -6.000  16.200


######################################################
def read_pdbid( input_pdb_file ):

	this_pdbid = input_pdb_file[0:4]
	this_pdbid = this_pdbid.lower()
	return this_pdbid


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
def read_PDB_file( pdb_file ):

	global got_DUM_N, got_DUM_O, DUM_N_z, DUM_O_z

	infile = open( pdb_file, "r" )
	inlines = infile.readlines()
	infile.close()

	prev_chain = ''
	prev_resn = ''
	prev_resi = 0
	prev_atom = ''
	got_DUM_N = 0
	got_DUM_O = 0
	DUM_N_z = 0
	DUM_O_z = 0

	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			if (len(inline) >= 6):
				this_key = inline[0:6]
				if ((this_key == 'HETATM') or (this_key == 'ATOM  ')):
					this_resn = inline[17:20]
					this_chain = inline[21:22]
					this_atom = inline[13:16]
					this_atom = this_atom.strip()
					this_resi = inline[22:26]
					this_resi = int(this_resi)
					this_z = inline[46:54]
					this_z = float(this_z)
					if (this_resn == 'DUM'):
						if ((got_DUM_N == 0) and (this_atom == 'N')):
							DUM_N_z = this_z
							got_DUM_N = 1
						if ((got_DUM_O == 0) and (this_atom == 'O')):
							DUM_O_z = this_z
							got_DUM_O = 1
					else:
						if (this_key == 'ATOM  '):
							this_is_residue = is_residue( this_resn )
							if (this_is_residue == 1):
								#same_as_prev_residue = 0
								#if ((this_chain == prev_chain) and (this_resi == prev_resi)):
								#	same_as_prev_residue = 1
								#if (same_as_prev_residue == 0):
								if (this_atom == 'CA'):
									this_chain_resi_already_exists = 0
									for i in range( 0, len(chain_array) ):
										compare_chain = chain_array[i]
										compare_resi = resi_array[i]
										if ((compare_chain == this_chain) and (compare_resi == this_resi)):
											this_chain_resi_already_exists = 1
									if (this_chain_resi_already_exists == 0):
										chain_array.append( this_chain )
										resi_array.append( this_resi )
										resn_array.append( this_resn )
										z_array.append( this_z )
										prev_chain = this_chain
										prev_resi = this_resi
	return


######################################################
def assign_positions():

	global got_DUM_N, got_DUM_O, DUM_N_z, DUM_O_z

	for i in range( 0, len(chain_array) ):
		this_chain = chain_array[i]
		this_resi = resi_array[i]
		this_resn = resn_array[i]
		this_z = z_array[i]
		posn = 'U' # unknown
		if ((got_DUM_N == 1) and (got_DUM_O == 1)):
			if (DUM_O_z > DUM_N_z):
				if ((this_z <= DUM_O_z) and (this_z >= DUM_N_z)):
					posn = 'M' # membrane
				elif (this_z >= DUM_O_z):
					posn = 'O' # outside (on the cytoplasmic side (or mitochondrial intermembrane space))
				else: # (this_z <= DUM_N_z)
					posn = 'I' # inside (on the periplasmic side (or mitochondrial intermembrane space))
			else: # (DUM_N_z > DUM_O_z):
				if ((this_z <= DUM_N_z) and (this_z >= DUM_O_z)):
					posn = 'M' # membrane
				elif (this_z >= DUM_N_z):
					posn = 'I' # inside (on the cytoplasmic side (or mitochondrial intermembrane space))
				else: # (this_z <= DUM_O_z)
					posn = 'O' # outside (on the periplasmic side (or mitochondrial intermembrane space))
		elif ((got_DUM_N == 0) and (got_DUM_O == 1)):
			if (this_z >= DUM_O_z):
				posn = 'O' # outside (on the periplasmic side (or mitochondrial intermembrane space))
			else: # (this_z < DUM_O_z)
				posn = 'M' # membrane
		elif ((got_DUM_N == 1) and (got_DUM_O == 0)):
			if (this_z <= DUM_N_z):
				posn = 'I' # inside (on the cytoplasmic side (or mitochondrial intermembrane space))
			else: # (this_z > DUM_N_z)
				posn = 'M' # membrane
		posn_array.append( posn )

		if (posn == 'M'):
			distance_to_outside = abs(DUM_O_z - this_z)
			distance_to_inside = abs(this_z - DUM_N_z)
			if (distance_to_outside >= distance_to_inside):
				membrane_near = 'I'
			else:
				membrane_near = 'O'
			membrane_near_array.append( membrane_near )
		else:
			membrane_near_array.append( '_' )

	#print_chain_array = ''
	#print_posn_array = ''
	#print_membrane_near_array = ''
	#for i in range( 0, len(chain_array) ):
	#	print_chain_array = print_chain_array + chain_array[i]
	#for i in range( 0, len(posn_array) ):
	#	print_posn_array = print_posn_array + posn_array[i]
	#for i in range( 0, len(membrane_near_array) ):
	#	print_membrane_near_array = print_membrane_near_array + membrane_near_array[i]
	#print print_chain_array
	#print print_posn_array
	#print print_membrane_near_array

	return


######################################################
def derive_segments_from_positions():

	prev_chain = ''
	prev_posn = ''
	prev_resi = 0
	concat_membrane_segments = ''
	concat_inside_segments = ''
	concat_outside_segments = ''
	open_membrane_segment = 0
	open_inside_segment = 0
	open_outside_segment = 0

	for i in range( 0, len(chain_array) ):
		this_chain = chain_array[i]
		this_resi = resi_array[i]
		this_resi = int(this_resi)
		this_resi_minus_one = this_resi - 1
		this_resn = resn_array[i]
		this_z = z_array[i]
		this_posn = posn_array[i]

		if (this_chain != prev_chain):

			# close off the segments for the previous chain (if there is one)
			if (prev_chain != ''):
				this_to = prev_resi
				if (open_membrane_segment == 1):
					concat_membrane_segments = concat_membrane_segments + str(this_to) + ','
					open_membrane_segment = 0
				elif (open_inside_segment == 1):
					concat_inside_segments = concat_inside_segments + str(this_to) + ','
					open_inside_segment = 0
				elif (open_outside_segment == 1):
					concat_outside_segments = concat_outside_segments + str(this_to) + ','
					open_outside_segment = 0
				membrane_segments_chain.append( prev_chain )
				membrane_segments_segments.append( concat_membrane_segments )
				inside_segments_chain.append( prev_chain )
				inside_segments_segments.append( concat_inside_segments )
				outside_segments_chain.append( prev_chain )
				outside_segments_segments.append( concat_outside_segments )

				prev_chain = this_chain
				prev_posn = ''
				concat_membrane_segments = ''
				concat_inside_segments = ''
				concat_outside_segments = ''

			# open the new segment of the new chain
			if (this_posn == 'M'):
				open_membrane_segment = 1
				concat_membrane_segments = str(this_resi) + '-'
			elif (this_posn == 'I'):
				open_inside_segment = 1
				concat_inside_segments = str(this_resi) + '-'
			elif (this_posn == 'O'):
				open_outside_segment = 1
				concat_outside_segments = str(this_resi) + '-'

 
		elif ((this_posn != prev_posn) or (prev_resi != this_resi_minus_one)):

			# close prev segment (if there is one)
			if (prev_posn != ''):
				this_to = prev_resi
				if (open_membrane_segment == 1):
					concat_membrane_segments = concat_membrane_segments + str(this_to) + ','
					open_membrane_segment = 0
				elif (open_inside_segment == 1):
					concat_inside_segments = concat_inside_segments + str(this_to) + ','
					open_inside_segment = 0
				elif (open_outside_segment == 1):
					concat_outside_segments = concat_outside_segments + str(this_to) + ','
					open_outside_segment = 0

			# open new segment in current chain
			if (this_posn == 'M'):
				open_membrane_segment = 1
				concat_membrane_segments = concat_membrane_segments + str(this_resi) + '-'
			elif (this_posn == 'I'):
				open_inside_segment = 1
				concat_inside_segments = concat_inside_segments + str(this_resi) + '-'
			elif (this_posn == 'O'):
				open_outside_segment = 1
				concat_outside_segments = concat_outside_segments + str(this_resi) + '-'

		prev_chain = this_chain
		prev_posn = this_posn
		prev_resi = this_resi

	# close off the segments for the last chain
	this_to = prev_resi
	if (open_membrane_segment == 1):
		concat_membrane_segments = concat_membrane_segments + str(this_to) + ','
	elif (open_inside_segment == 1):
		concat_inside_segments = concat_inside_segments + str(this_to) + ','
	elif (open_outside_segment == 1):
		concat_outside_segments = concat_outside_segments + str(this_to) + ','
	membrane_segments_chain.append( prev_chain )
	membrane_segments_segments.append( concat_membrane_segments )
	inside_segments_chain.append( prev_chain )
	inside_segments_segments.append( concat_inside_segments )
	outside_segments_chain.append( prev_chain )
	outside_segments_segments.append( concat_outside_segments )

	return


######################################################
def collate_segments_by_chain():

	for this_chain in membrane_segments_chain:
		chain_exists = 0
		for existing_chain in all_segments_chain:
			if (this_chain == existing_chain):
				chain_exists = 1
		if (chain_exists == 0):
			all_segments_chain.append( this_chain )
			all_segments_segments.append( '' )
	for this_chain in inside_segments_chain:
		chain_exists = 0
		for existing_chain in all_segments_chain:
			if (this_chain == existing_chain):
				chain_exists = 1
		if (chain_exists == 0):
			all_segments_chain.append( this_chain )
			all_segments_segments.append( '' )
	for this_chain in outside_segments_chain:
		chain_exists = 0
		for existing_chain in all_segments_chain:
			if (this_chain == existing_chain):
				chain_exists = 1
		if (chain_exists == 0):
			all_segments_chain.append( this_chain )
			all_segments_segments.append( '' )

	for i in range( 0, len(membrane_segments_segments) ):
		this_segments = membrane_segments_segments[i]
		this_chain = membrane_segments_chain[i]
		if (this_segments != ''):
			for j in range( 0, len(all_segments_chain) ):
				existing_chain = all_segments_chain[j]
				if (existing_chain == this_chain):
					existing_segments = all_segments_segments[j]
					existing_segments = existing_segments + 'MEMBRANE-SEGMENTS,,,,,' + this_segments + ';;;;;'
					all_segments_segments[j] = existing_segments
	for i in range( 0, len(inside_segments_segments) ):
		this_segments = inside_segments_segments[i]
		this_chain = inside_segments_chain[i]
		if (this_segments != ''):
			for j in range( 0, len(all_segments_chain) ):
				existing_chain = all_segments_chain[j]
				if (existing_chain == this_chain):
					existing_segments = all_segments_segments[j]
					existing_segments = existing_segments + 'INSIDE-SEGMENTS,,,,,' + this_segments + ';;;;;'
					all_segments_segments[j] = existing_segments
	for i in range( 0, len(outside_segments_segments) ):
		this_segments = outside_segments_segments[i]
		this_chain = outside_segments_chain[i]
		if (this_segments != ''):
			for j in range( 0, len(all_segments_chain) ):
				existing_chain = all_segments_chain[j]
				if (existing_chain == this_chain):
					existing_segments = all_segments_segments[j]
					existing_segments = existing_segments + 'OUTSIDE-SEGMENTS,,,,,' + this_segments + ';;;;;'
					all_segments_segments[j] = existing_segments

	return


######################################################
def write_output_segments( this_pdbid, output_file_name ):

	outfile = open( output_file_name, "a" ) # open output file for writing in append mode

	for i in range( 0, len(all_segments_chain) ):
		this_chain = all_segments_chain[i]
		this_segments = all_segments_segments[i]
		this_id = this_pdbid + '_' + this_chain
		output_line = this_id + ':::::' + this_segments + ':::::'
		output_line = output_line + "\r\n"
		outfile.write( output_line )

	outfile.close()

	return


######################################################
def look_for_item( an_item, an_array ):

	see_this_item = 0
	for this_item in an_array:
		if (this_item == an_item):
			see_this_item = 1

	return see_this_item


######################################################
def write_output_seq( this_pdbid, output_file_name ):

	chainseq_chain_array = []
	chainseq_seq_array = []

	prev_chain = ''
	prev_resi = 0
	prev_resn = ''
	seq_string = ''

	for i in range( 0, len(chain_array) ):
		this_chain = chain_array[i]
		this_resi = resi_array[i]
		this_resi = int(this_resi)
		this_resn = resn_array[i]

		if (this_chain != prev_chain):

			# at the end of reading a chain, save the entire chain sequence
			if (prev_chain != ''):
				seen_this_chain_before = look_for_item( prev_chain, chainseq_chain_array )
				# if we've already seen this chain before, then this second appearance is probably an amino acid ligand entered as 'ATOM  '
				if (seen_this_chain_before == 0):
					chainseq_chain_array.append( prev_chain )
					chainseq_seq_array.append( seq_string )

			prev_resi = 0
			seq_string = ''
			#if (this_resi > 1):
			#	for j in range( 0, (this_resi - 1) ):
			#		seq_string = seq_string + '_'

		if ((prev_resi + 1) < this_resi):
			for j in range( (prev_resi + 1), this_resi ):
				seq_string = seq_string + '_'
		this_resn_1char = resn_3char_to_1char( this_resn )
		seq_string = seq_string + this_resn_1char
		prev_chain = this_chain
		prev_resi = this_resi

	# at the end of reading a chain, save the entire chain sequence
	if (prev_chain != ''):
		seen_this_chain_before = look_for_item( prev_chain, chainseq_chain_array )
		# if we've already seen this chain before, then this second appearance is probably an amino acid ligand entered as 'ATOM  '
		if (seen_this_chain_before == 0):
			chainseq_chain_array.append( prev_chain )
			chainseq_seq_array.append( seq_string )

	# now that we have nice arrays for chains and sequences, write to output file

	outfile = open( output_file_name, "a" ) # open output file for writing in append mode

	for i in range( 0, len(chainseq_chain_array) ):
		this_chain = chainseq_chain_array[i]
		this_sequence = chainseq_seq_array[i]
		this_id = this_pdbid + '_' + this_chain
		output_line = this_id + ':::::' + this_sequence + ':::::'
		output_line = output_line + "\r\n"
		outfile.write( output_line )

	outfile.close()

	return


######################################################
def write_output_membrane_near( this_pdbid, output_file_name ):

	chainnear_chain_array = []
	chainnear_membrane_near_array = []

	prev_chain = ''
	prev_resi = 0
	prev_membrane_near = ''
	membrane_near_string = ''

	for i in range( 0, len(chain_array) ):
		this_chain = chain_array[i]
		this_resi = resi_array[i]
		this_membrane_near = membrane_near_array[i]

		if (this_chain != prev_chain):

			# at the end of reading a chain, save the entire chain sequence
			if (prev_chain != ''):
				seen_this_chain_before = look_for_item( prev_chain, chainnear_chain_array )
				# if we've already seen this chain before, then this second appearance is probably an amino acid ligand entered as 'ATOM  '
				if (seen_this_chain_before == 0):
					chainnear_chain_array.append( prev_chain )
					chainnear_membrane_near_array.append( membrane_near_string )

			prev_resi = 0
			prev_membrane_near = ''
			membrane_near_string = ''

		if ((prev_resi + 1) < this_resi):
			for j in range( (prev_resi + 1), this_resi ):
				membrane_near_string = membrane_near_string + '_'
		membrane_near_string = membrane_near_string + this_membrane_near
		prev_chain = this_chain
		prev_resi = this_resi
		prev_membrane_near = this_membrane_near

	# at the end of reading a chain, save the entire chain sequence
	if (prev_chain != ''):
		seen_this_chain_before = look_for_item( prev_chain, chainnear_chain_array )
		# if we've already seen this chain before, then this second appearance is probably an amino acid ligand entered as 'ATOM  '
		if (seen_this_chain_before == 0):
			chainnear_chain_array.append( prev_chain )
			chainnear_membrane_near_array.append( membrane_near_string )

	# now that we have nice arrays for chains and sequences, write to output file

	outfile = open( output_file_name, "a" ) # open output file for writing in append mode

	for i in range( 0, len(chainnear_chain_array) ):
		this_chain = chainnear_chain_array[i]
		this_membrane_near = chainnear_membrane_near_array[i]
		this_id = this_pdbid + '_' + this_chain
		output_line = this_id + ':::::' + this_membrane_near + ':::::'
		output_line = output_line + "\r\n"
		outfile.write( output_line )

	outfile.close()

	return


######################################################
def main():

	global outfile
	input_pdb_file = sys.argv[1]
	output_file_name_segments = sys.argv[2]
	output_file_name_seq = sys.argv[3]
	output_file_name_membrane_near = sys.argv[4]

	this_pdbid = read_pdbid( input_pdb_file )

	read_PDB_file( input_pdb_file )

	assign_positions()

	#print 'chain_array =', chain_array
	#print 'resi_array =', resi_array
	#print 'resn_array =', resn_array
	#print 'z_array =', z_array
	#print 'posn_array =', posn_array
	#print 'membrane_near_array =', membrane_near_array

	derive_segments_from_positions()

	#print 'membrane_segments_segments =', membrane_segments_segments
	#print 'membrane_segments_chain = ', membrane_segments_chain
	#print 'inside_segments_segments = ', inside_segments_segments
	#print 'inside_segments_chain = ', inside_segments_chain
	#print 'outside_segments_segments = ', outside_segments_segments
	#print 'outside_segments_chain = ', outside_segments_chain

	collate_segments_by_chain()

	#print 'all_segments_segments = ', all_segments_segments
	#print 'all_segments_chain = ', all_segments_chain

	write_output_segments( this_pdbid, output_file_name_segments )

	write_output_seq( this_pdbid, output_file_name_seq )

	write_output_membrane_near( this_pdbid, output_file_name_membrane_near )

	# later, find coils, half-membrane helices and hairpins


if __name__ == "__main__":
	# Someone is launching this directly
	main()

